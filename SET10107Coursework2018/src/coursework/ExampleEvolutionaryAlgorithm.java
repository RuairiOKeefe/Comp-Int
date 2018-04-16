package coursework;

import java.util.ArrayList;

import model.Fitness;
import model.Individual;
import model.LunarParameters.DataSet;
import model.NeuralNetwork;

/**
 * Implements a basic Evolutionary Algorithm to train a Neural Network
 * 
 * You Can Use This Class to implement your EA or implement your own class that
 * extends {@link NeuralNetwork}
 * 
 */
public class ExampleEvolutionaryAlgorithm extends NeuralNetwork
{

	int firstParent = -1;

	/**
	 * The Main Evolutionary Loop
	 */
	@Override
	public void run()
	{

		// Initialise a population of Individuals with random weights
		population = initialise();

		// Record a copy of the best Individual in the population
		best = getBest();
		System.out.println("Best From Initialisation " + best);

		/**
		 * main EA processing loop
		 */

		while (evaluations < Parameters.maxEvaluations)
		{

			/**
			 * this is a skeleton EA - you need to add the methods. You can also change the
			 * EA if you want You must set the best Individual at the end of a run
			 * 
			 */

			// Select 2 Individuals from the current population.
			Individual parent1 = select();
			Individual parent2 = select();

			// Generate a child by crossover. Not Implemented
			ArrayList<Individual> children = reproduce(parent1, parent2);

			// mutate the offspring
			mutate(children);

			// Evaluate the children
			evaluateIndividuals(children);

			// Replace children in population
			replace(children);

			// check to see if the best has improved
			best = getBest();

			// Implemented in NN class.
			outputStats();

			// Increment number of completed generations
		}

		// save the trained network to disk
		saveNeuralNetwork();
	}

	/**
	 * Sets the fitness of the individuals passed as parameters (whole population)
	 * 
	 */
	private void evaluateIndividuals(ArrayList<Individual> individuals)
	{
		for (Individual individual : individuals)
		{
			individual.fitness = Fitness.evaluate(individual, this);
		}
	}

	/**
	 * Returns a copy of the best individual in the population
	 * 
	 */
	private Individual getBest()
	{
		best = null;
		;
		for (Individual individual : population)
		{
			if (best == null)
			{
				best = individual.copy();
			}
			else if (individual.fitness < best.fitness)
			{
				best = individual.copy();
			}
		}
		return best;
	}

	/**
	 * Generates a randomly initialised population
	 * 
	 */
	private ArrayList<Individual> initialise()
	{
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.popSize; ++i)
		{
			// chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			population.add(individual);
		}
		evaluateIndividuals(population);
		return population;
	}

	private int tournament(int first, int second)
	{
		if (population.get(first).compareTo(population.get(second)) == 1)
			return second;
		else
			return first;
	}

	private Individual tournamentSelect()
	{
		ArrayList<Integer> indexArray = new ArrayList<Integer>();
		for (int i = 0; i < population.size(); i++)
		{
			indexArray.add(i);
		}

		// Keep parents unique
		if (firstParent != -1)
		{
			indexArray.remove(firstParent);
		}
		int winner = indexArray.get((int) (Parameters.random.nextDouble() * indexArray.size()));

		indexArray.remove((int) (Parameters.random.nextDouble() * indexArray.size()));

		for (int i = 0; i < Parameters.tSize; i++)
		{

			int randomId = indexArray.get((int) (Parameters.random.nextDouble() * indexArray.size()));

			indexArray.remove((int) (Parameters.random.nextDouble() * indexArray.size()));

			winner = tournament(winner, randomId);

		}
		Individual parent = population.get(winner);
		firstParent = winner;

		return parent.copy();
	}

	private Individual select()
	{

		Individual parent = population.get(Parameters.random.nextInt(Parameters.popSize));
		switch (Parameters.selectType)
		{
		case RANDOM:
			break;
		case TOURNAMENT:
			parent = tournamentSelect();
			break;
		default:
			break;
		}
		return parent.copy();
	}

	private ArrayList<Individual> uniformCross(Individual parent1, Individual parent2)
	{
		ArrayList<Individual> children = new ArrayList<>();
		Individual child = new Individual();
		double[] chromosome = new double[Parameters.getNumGenes()];
		for (int i = 0; i < Parameters.getNumGenes(); i++)
		{
			if (Parameters.random.nextDouble() < Parameters.uniformProb)
			{
				chromosome[i] = parent1.chromosome[i];
			}
			else
			{
				chromosome[i] = parent2.chromosome[i];
			}
		}
		child.chromosome = chromosome;
		children.add(child);
		return children;
	}

	private ArrayList<Individual> onePointCross(Individual parent1, Individual parent2)
	{
		ArrayList<Individual> children = new ArrayList<>();
		Individual child = new Individual();
		boolean parent = false;
		double[] chromosome = new double[Parameters.getNumGenes()];
		int flipPoint = (int) (Parameters.random.nextDouble() * Parameters.getNumGenes());
		for (int i = 0; i < Parameters.getNumGenes(); i++)
		{
			if (i > flipPoint)
				parent = !parent;

			if (parent)
			{
				chromosome[i] = parent1.chromosome[i];
			}
			else
			{
				chromosome[i] = parent2.chromosome[i];
			}
		}
		child.chromosome = chromosome;
		children.add(child);
		return children;
	}

	private ArrayList<Individual> nPointCross(Individual parent1, Individual parent2)
	{
		ArrayList<Individual> children = new ArrayList<>();
		Individual child = new Individual();
		boolean parent = false;
		double[] chromosome = new double[Parameters.getNumGenes()];
		ArrayList<Integer> flipPoint = new ArrayList<Integer>();

		for (int i = 0; i < Parameters.numPoints; i++)
		{
			flipPoint.add((int) (Parameters.random.nextDouble() * Parameters.getNumGenes()));
		}

		for (int i = 0; i < Parameters.getNumGenes(); i++)
		{
			for (int j = 0; j < flipPoint.size(); j++)
			{
				if (i > flipPoint.get(j))
				{
					parent = !parent;
					flipPoint.remove(j);
				}
			}
			if (parent)
			{
				chromosome[i] = parent1.chromosome[i];
			}
			else
			{
				chromosome[i] = parent2.chromosome[i];
			}
		}
		child.chromosome = chromosome;
		children.add(child);
		return children;
	}

	private ArrayList<Individual> reproduce(Individual parent1, Individual parent2)
	{
		ArrayList<Individual> children = new ArrayList<>();
		switch (Parameters.crossoverType)
		{
		case CLONE:
			children.add(parent1.copy());
			children.add(parent2.copy());
			break;
		case UNIFORM:
			children = uniformCross(parent1, parent2);
			break;
		case ONEPOINT:
			children = onePointCross(parent1, parent2);
			break;
		case NPOINT:
			children = nPointCross(parent1, parent2);
			break;
		}
		return children;
	}

	/**
	 * Mutation
	 * 
	 * 
	 */
	private void mutate(ArrayList<Individual> individuals)
	{
		for (Individual individual : individuals)
		{
			for (int i = 0; i < individual.chromosome.length; i++)
			{
				if (Parameters.random.nextDouble() < Parameters.mutateRate)
				{
					if (Parameters.random.nextBoolean())
					{
						individual.chromosome[i] += (Parameters.mutateChange);
					}
					else
					{
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
		}
	}

	private void randomReplace(ArrayList<Individual> individuals)
	{
		for (int i = 0; i < Parameters.numReplacements; i++)
		{
			int randParent = (int) (Parameters.random.nextDouble() * population.size());
			int randChild = (int) (Parameters.random.nextDouble() * individuals.size());
			population.set(randParent, individuals.get(randChild));
		}
	}

	private void tournamentReplace(ArrayList<Individual> individuals)
	{
		ArrayList<Individual> newPop = new ArrayList<Individual>();
		for (int n = 0; n < Parameters.popSize; n++)
		{
			ArrayList<Individual> combinedPop = population;

			for (int i = 0; i < individuals.size(); i++)
			{
				combinedPop.add(individuals.get(i));
			}
			ArrayList<Integer> indexArray = new ArrayList<Integer>();
			for (int i = 0; i < combinedPop.size(); i++)
			{
				indexArray.add(i);
			}

			int winner = indexArray.get((int) (Parameters.random.nextDouble() * indexArray.size()));

			indexArray.remove((int) (Parameters.random.nextDouble() * indexArray.size()));

			for (int i = 0; i < Parameters.rtSize; i++)
			{

				int randomId = indexArray.get((int) (Parameters.random.nextDouble() * indexArray.size()));

				indexArray.remove((int) (Parameters.random.nextDouble() * indexArray.size()));

				winner = tournament(winner, randomId);

				if (combinedPop.get(winner).compareTo(combinedPop.get(randomId)) != 1)
					winner = randomId;
			}
			Individual newIndividual = combinedPop.get(winner);
			newPop.add(newIndividual);
		}
		population = newPop;
	}

	private void replaceWorst(ArrayList<Individual> individuals)
	{
		for (Individual individual : individuals)
		{
			int idx = getWorstIndex();
			population.set(idx, individual);
		}
	}

	private void replace(ArrayList<Individual> individuals)
	{
		switch (Parameters.replaceType)
		{
		case RANDOM:
			randomReplace(individuals);
			break;
		case TOURNAMENT:
			tournamentReplace(individuals);
			break;
		case WORST:
			replaceWorst(individuals);
			break;
		}
	}

	/**
	 * Returns the index of the worst member of the population
	 * 
	 * @return
	 */
	private int getWorstIndex()
	{
		Individual worst = null;
		int idx = -1;
		for (int i = 0; i < population.size(); i++)
		{
			Individual individual = population.get(i);
			if (worst == null)
			{
				worst = individual;
				idx = i;
			}
			else if (individual.fitness > worst.fitness)
			{
				worst = individual;
				idx = i;
			}
		}
		return idx;
	}

	@Override
	public double activationFunction(double x)
	{
		if (x < -20.0)
		{
			return -1.0;
		}
		else if (x > 20.0)
		{
			return 1.0;
		}
		return Math.tanh(x);
	}
}
