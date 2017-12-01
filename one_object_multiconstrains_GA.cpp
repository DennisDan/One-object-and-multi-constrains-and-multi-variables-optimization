#include <bitset>
#include<iostream>
#include<string>
#include<vector>
#include<time.h>
#include<sstream>
using namespace std;



#define population_number 100
#define gene_length 80  // 40 -40   the longer length the more precise
#define child_number 100
#define iteration  100
#define crossover_rate  0.5   //crossover rate and mutation rate decide the results
#define mutation_rate  0.25
#define penalty_coefficient 10



#define GAposrand() (rand()%(gene_length-1))  //0~2  generate a integer
#define BinaryRand() (rand()%2)  //generate 1 or 0 
#define Srand() (float)rand()/(float)RAND_MAX  //0~1 generate a float


// boundary conditions
#define x1_up  100   
#define x1_low 13
#define x2_up  100 
#define x2_low  0





//  binary transformation
long double binary_transform(string str)
{
	int num = 0;
	int val;
	string itr;
	for (size_t i = 0; i < str.length(); i++)
	{
		itr = str[i];
		istringstream iss(itr);
		iss >> val;
		num += val*pow(2, str.length() - 1 - i);
	}

	return num;
}




int Binary2Hex(std::string Binary)
{
	std::bitset<32> set(Binary);
	int hex = set.to_ulong();

	return hex;
}

// Convert the 32-bit binary into the decimal
float GetFloat32(std::string Binary)
{
	int HexNumber = Binary2Hex(Binary);

	bool negative = !!(HexNumber & 0x80000000);
	int  exponent = (HexNumber & 0x7f800000) >> 23;
	int sign = negative ? -1 : 1;

	// Subtract 127 from the exponent
	exponent -= 127;

	// Convert the mantissa into decimal using the
	// last 23 bits
	int power = -1;
	float total = 0.0;
	for (int i = 0; i < 23; i++)
	{
		int c = Binary[i + 9] - '0';
		total += (float)c * (float)pow(2.0, power);
		power--;
	}
	total += 1.0;

	float value = sign * (float)pow(2.0, exponent) * total;

	return value;
}

// Get 32-bit IEEE 754 format of the decimal value
std::string GetBinary32(float value)
{
	union
	{
		float input;   // assumes sizeof(float) == sizeof(int)
		int   output;
	}    data;

	data.input = value;

	std::bitset<sizeof(float) * CHAR_BIT>   bits(data.output);

	std::string mystring = bits.to_string<char,
		std::char_traits<char>,
		std::allocator<char> >();

	return mystring;
}



typedef struct PARENT
{
	long double value;  //for fitness  倒數過了
	int genes[gene_length];
	long double fitness;
	long double x1, x2;
	long double f;


}parent;




class GA {

	parent population[population_number];
	parent child[child_number];
	parent total[population_number + child_number];
public:

	void initialize();
	void crossover();   //用 do while試試
	void mutation();
	void selection();
	void fitness();
	void binary_conversion();
	void print();
};

void GA::print() {

	cout << " value  = " << total[0].f << "  x1  x2   " << total[0].x1 << "    " << total[0].x2 << endl;


}

void GA::initialize() {

	for (size_t i = 0; i < population_number; i++)
	{
		for (size_t j = 0; j < gene_length; j++)
		{
			population[i].genes[j] = BinaryRand();
		}
	}
}

void GA::crossover() {

	for (size_t i = 0; i < child_number / 2; i++)
	{
		int crossPoint = GAposrand();
		int parent_a = (rand() % (population_number - 1));
		int parent_b = (rand() % (population_number - 1));



		if (Srand() > crossover_rate)
		{

			for (size_t j = 0; j < crossPoint; j++)
			{
				child[i * 2].genes[j] = population[parent_a].genes[j];
				child[i * 2 + 1].genes[j] = population[parent_b].genes[j];

			}

			for (size_t k = crossPoint; k < gene_length; k++)
			{
				child[i * 2].genes[k] = population[parent_b].genes[k];
				child[i * 2 + 1].genes[k] = population[parent_a].genes[k];
			}
		}

		else
		{
			child[i * 2] = population[parent_a];
			child[i * 2 + 1] = population[parent_b];

		}


	}
}

void GA::mutation()
{
	for (size_t i = 0; i < child_number; i++)
	{
		for (size_t j = 0; j < gene_length; j++)
		{

			if (Srand() < mutation_rate)
			{
				if (child[i].genes[j] == 0)child[i].genes[j] = 1;
				else child[i].genes[j] = 0;
			}
		}
	}
}

void GA::binary_conversion() {


	for (size_t i = 0; i < population_number; i++)
	{
		memcpy(&total[i], &population[i], sizeof(parent));
	}

	for (size_t i = population_number; i < population_number + child_number; i++)
	{
		memcpy(&total[i], &child[i], sizeof(parent));
	}

	for (size_t i = 0; i < population_number + child_number; i++)
	{
		string aaa, bbb;

		for (size_t j = 0; j < gene_length / 2; j++)    
		{
			aaa += std::to_string(total[i].genes[j]);
		}
		for (size_t k = gene_length / 2; k < gene_length; k++)    
		{
			bbb += std::to_string(total[i].genes[k]);
		}

		total[i].x1 = x1_low + binary_transform(aaa)*(x1_up - x1_low) / (pow(2, gene_length / 2) - 1);
		total[i].x2 = x2_low + binary_transform(bbb)*(x2_up - x2_low) / (pow(2, gene_length / 2) - 1);

	}
}

void GA::fitness() {


	for (size_t i = 0; i < population_number + child_number; i++)
	{

		if ((-pow((total[i].x1 - 5), 2) - pow((total[i].x2 - 5), 2) + 100) <= 0 &&
			(pow(total[i].x1 - 6, 2) + pow(total[i].x2 - 5, 2) - 82.81) <= 0
			)
		{
			total[i].value = 1 / (pow((total[i].x1 - 10), 3) + pow((total[i].x2 - 20), 3));

		}

		else
		{
			total[i].value =
				1 / (
					pow((total[i].x1 - 10), 3) + pow((total[i].x2 - 20), 3) +
					1 * penalty_coefficient*
					(
						pow((-pow((total[i].x1 - 5), 2) - pow((total[i].x2 - 5), 2) + 100), 2) +
						pow((pow(total[i].x1 - 6, 2) + pow(total[i].x2 - 5, 2) - 82.81), 2)

						)
					);

		}
		total[i].f = pow((total[i].x1 - 10), 3) + pow(total[i].x2 - 20, 3);

		cout << " value  = " << total[i].f << "  x1  x2   " << total[i].x1 << "    " << total[i].x2 << endl;
	}
};

void GA::selection()
{
	long double total_fitness = 0;
	double probability[population_number + child_number];


	for (size_t i = 0; i < population_number + child_number; i++)
	{
		total_fitness += total[i].value;
	}

	for (size_t j = 0; j < population_number + child_number; j++)
	{
		probability[j] = total[j].value / total_fitness;
	}


	for (size_t k = 0; k < population_number; k++)
	{
		double criteria = Srand();
		double prob = 0;
		int w = 0;
		for (size_t w = 0; w < population_number + child_number; w++)
		{
			prob += probability[w];
			if (prob >= criteria)
			{
				break;
			}


		}

		memcpy(&population[k], &total[w], sizeof(parent));
	}
}


int main()
{
	srand((unsigned)time(NULL));  
	GA asd;
	asd.initialize();

	for (size_t i = 0; i < iteration; i++)
	{
		asd.crossover();
		asd.mutation();
		asd.binary_conversion();
		asd.fitness();
		asd.selection();
	}



	system("pause");
	return 0;
}


















