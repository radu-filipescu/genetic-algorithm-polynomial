#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const int MUTATION_MODE = 0; // 0 - rare mutation, 1 - common mutation

int populationSize;
double domainLeft, domainRight;
double a, b, c; // polynomial coefficients
double precision;
double crossoverChance;
double mutationChance;
int generations;
mt19937 rng;
ofstream fout("results.txt");

class Individual {
public:
    vector <int> chromosome;
    double fitness;

    Individual() {
        chromosome.resize(precision);
    }
};

vector <Individual> population;

void readFromInput() {
    ifstream fin("input.txt");

    fin >> populationSize;
    fin >> domainLeft >> domainRight;
    fin >> a >> b >> c;
    fin >> precision;
    fin >> crossoverChance;
    fin >> mutationChance;
    fin >> generations;

    fin.close();
}

// wrapper for <random> library mt19937
// returns a random number [0,1)
double randomGenerator() {
    return (double)rng() / ((long long)rng.max() + 1);
}

double getValueFromChromosome(vector <int> chromosome) {
    double x = ( (domainRight - domainLeft) / (pow(2, chromosome.size()) - 1));
    double ans = 0;
    for(int i = chromosome.size() - 1; i >= 0; --i) {
        ans += chromosome[i] * x;
        x *= 2;
    }
    return ans + domainLeft;
}

void computeFitness() {
    for(int i = 0; i < population.size(); ++i) {
        double variableValue = getValueFromChromosome(population[i].chromosome);
        population[i].fitness = a * variableValue * variableValue + b * variableValue + c;
    }
}

void initialize() {
    rng.seed(chrono::steady_clock::now().time_since_epoch().count());

    long long chromosomeLength = 1;
    for(int i = 0; i < precision; i++)
        chromosomeLength *= 10;
    chromosomeLength *= (domainRight - domainLeft);
    chromosomeLength = log2(chromosomeLength);

    // create population
    for(int i = 0; i < populationSize; i++){
        Individual tmp;
        tmp.chromosome.resize(chromosomeLength);

        // initial generation has random genes
        for(int j = 0; j < tmp.chromosome.size(); ++j)
            tmp.chromosome[j] = randomGenerator() / (double)0.5;

        population.push_back(tmp);
    }

    // compute initial fitness
    computeFitness();
}

int binSearch(int lf, int rg, double val, vector<double> &v) {
    if(lf > rg) return -1;

    int mid = (lf + rg) / 2;

    if(v[mid] < val)
        return max(mid, binSearch(mid + 1, rg, val, v));

    return binSearch(lf, mid - 1, val, v);
}

void selection() {
    vector <double> probs;

    double totalFitness = 0;
    for(int i = 0; i < population.size(); ++i)
        totalFitness += population[i].fitness;

    for(int i = 0; i < population.size(); ++i) {
        double currentProb = 0;
        if(i > 0) currentProb = probs[i - 1];

        probs.push_back(currentProb + population[i].fitness / totalFitness);
    }

    // generate random numbers and pick individuals
    vector <Individual> temporary;

    for(int i = 0; i < populationSize; ++i) {
        // pick an individual that goes in the next generation
        int pickedNumber = 1 + binSearch(0, probs.size() - 1, randomGenerator(), probs);

        temporary.push_back(population[pickedNumber]);
    }

    population = temporary;
}

void crossGenes(Individual & X, Individual & Y) {
    // pick a random point in the chromosome for crossover

    int pickedPoint = (randomGenerator() / 1) * X.chromosome.size();
    if(pickedPoint == 0)
        pickedPoint += 1;

    for(int i = 0; i < pickedPoint; ++i)
        swap(X.chromosome[i], Y.chromosome[i]);
}

void crossover() {
    // select which individuals will enter crossover stage ( based on crossover chance )
    vector <Individual> v;
    vector <bool> selected;

    for(int i = 0; i < population.size(); ++i) {
        selected.push_back(false);

        if(randomGenerator() <= crossoverChance) {
            v.push_back(population[i]);

            selected[i] = true;
        }
    }

    if(v.size() < 2) return;

    for(int i = 0; i < v.size() - 1; i += 2)
        crossGenes(v[i], v[i+1]);

    // replace the selected Individuals with their crossed-over versions
    int cursor = 0;
    for(int j = 0; cursor < v.size() && j < population.size(); ++j)
        if(selected[j]) {
            population[j] = v[cursor];
            ++cursor;
        }
}

void mutateIndividual(Individual & X) {
    // rare mutation
    if(MUTATION_MODE == 0) {
        if(randomGenerator() <= mutationChance) {
            int randomPosition = (randomGenerator() / 1) * X.chromosome.size();

            X.chromosome[randomPosition] = X.chromosome[randomPosition] ^ 1 ^ 0;
        }
    }
    // common mutation
    else {
        for(int i = 0; i < X.chromosome.size(); ++i)
            if(randomGenerator() <= mutationChance)
                X.chromosome[i] = X.chromosome[i] ^ 1 ^ 0;
    }
}

void mutation() {
    for(int i = 0; i < population.size(); ++i)
        mutateIndividual(population[i]);
}

// this function is made only for printing in the output file
// and illustrating the process, the resulting
// individuals won't be used in the actual simulation

void printInitialState() {
    vector <Individual> population2;

    population2 = population;

    fout << "Initial state:\n";

    for(int i = 0; i < population2.size(); ++i) {
        fout << "*" << ((i<10)?" " : "" )<< i << ":\n";
        fout << "chromosome:     ";
        for(int j = 0; j < population2[i].chromosome.size(); ++j)
            fout << population2[i].chromosome[j];
        fout << "\nvariable value: " << getValueFromChromosome(population2[i].chromosome);
        fout << "\nfitness:        " << population2[i].fitness << "\n\n";
    }

    fout << "selection chances:\n";

    vector <double> probs;

    double totalFitness = 0;
    for(int i = 0; i < population2.size(); ++i)
        totalFitness += population2[i].fitness;

    for(int i = 0; i < population2.size(); ++i) {
        double currentProb = 0;
        if(i > 0) currentProb = probs[i - 1];

        probs.push_back(currentProb + population2[i].fitness / totalFitness);

        fout << "*" << ((i<10)?" " : "" )<< i << ": ";
        fout << population2[i].fitness / totalFitness << '\n';
    }

    fout << "\nselection intervals:\n";
    for(int i = 0; i < population2.size(); ++i)
        fout << probs[i] << ' ';
    fout << "\n\n";

    // emulate a generation run for printing into file
    fout << "picking the next generation:\n";
    vector<Individual> temporary;
    for(int i = 0; i < populationSize; ++i) {
        double randomValue = randomGenerator();

        int pickedNumber = 1 + binSearch(0, probs.size() - 1, randomValue, probs);
        temporary.push_back(population2[pickedNumber]);

        fout << "random value: " << fixed << setprecision(5) << randomValue << ' ' << "picked number: " << pickedNumber << '\n';
    }
    population2 = temporary;

    fout << "\nnew generation:\n";
    for(int i = 0; i<population2.size(); ++i) {
        fout << "*" << ((i<10)?" " : "" )<< i << ": ";
        for(int j = 0; j <population2[i].chromosome.size(); ++j)
            fout << population2[i].chromosome[j];
        fout << '\n';
    }

    fout << "\ncrossover chance: " << crossoverChance << '\n';
    vector <Individual> v;
    vector <bool> selected;
    for(int i = 0; i < population2.size(); ++i) {
        selected.push_back(false);
        double randomValue = randomGenerator();
        if(randomValue <= crossoverChance) {
            v.push_back(population2[i]);

            selected[i] = true;
        }
        fout << "*" << ((i<10)?" " : "" )<< i << ": " << randomValue << ' ';

        if(selected[i]) fout << "gets selected";
        else fout << "doesn't get selected";

        fout << '\n';
    }
    fout << '\n';
    if(v.size() >= 2) {
        fout << "the swaps will be:\n";
        for(int i = 0; i < v.size() - 1; i += 2) {
            for(int j = 0; j < v[i].chromosome.size(); ++j)
                fout << v[i].chromosome[j];
            fout << " <-> ";
            for(int j = 0; j < v[i + 1].chromosome.size(); ++j)
                fout << v[i + 1].chromosome[j];
            fout << '\n';
            crossGenes(v[i], v[i+1]);
        }
        int cursor = 0;
        for(int j = 0; cursor < v.size() && j < population2.size(); ++j)
            if(selected[j]) {
                population2[j] = v[cursor];
                ++cursor;
            }
    }
    fout << "\nafter crossover:\n";
    for(int i = 0; i<population2.size(); ++i) {
        fout << "*" << ((i<10)?" " : "" )<< i << ": ";
        for(int j = 0; j <population2[i].chromosome.size(); ++j)
            fout << population2[i].chromosome[j];
        fout << '\n';
    }
    fout << "\nmutation chance: " << mutationChance << '\n';
    fout << "after mutations:\n";
    for(int i = 0; i < population2.size(); ++i) {
        fout << "*" << ((i<10)?" " : "" )<< i << ":\n";
        fout << "chromosome:     ";
        for(int j = 0; j < population2[i].chromosome.size(); ++j)
            fout << population2[i].chromosome[j];
        fout << "\nvariable value: " << getValueFromChromosome(population2[i].chromosome);
        fout << "\nfitness:        " << population2[i].fitness << "\n\n";
    }
}

int main()
{
    readFromInput();
    initialize();

    printInitialState();

    fout << "\nevolution:\n";

    for(int i = 0; i < generations; ++i) {
        selection();
        crossover();
        mutation();
        computeFitness();

        double maxFitness = -1;
        double maxValue;
        for(int j = 0; j < population.size(); ++j)
            if(population[j].fitness > maxFitness) {
                maxFitness = population[j].fitness;
                maxValue = getValueFromChromosome(population[j].chromosome);
            }
        fout << "value: " << maxFitness << " polynomial variable: " << maxValue << '\n';
    }

    return 0;
}
