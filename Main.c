#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/* + ----------------------------------------------------
   ! CONSTANTES DE LA FUNCION Y EL ALGORITMO GENÉTICO
   +-----------------------------------------------------*/
#define POPULATION_SIZE 100
#define FX_LOWER_BOUND -20
#define FX_UPPER_BOUND  20
#define PRECISION       3
/* + ----------------------------------------------------
   ! DEFINICION DEL INDIVIDUO
   +-----------------------------------------------------*/
typedef struct ind_t{
    int * chromosome;
    double x;                //Genes
    double fitness;          // El individuo va a terner una aptitude
    int parents[2];                         
    int mutation_place;
    int crossover_place;
}Individual;                // Ya esta la definicion para un solo individuo
/* + ----------------------------------------------------
   ! POBLACIONES (Y RULETA-->VIDEO2) (Y ELITISTA --VIDEO5)
   +-----------------------------------------------------*/
Individual* parents;        // Apuntador al arreglo de la población
Individual* offspring;
Individual the_best;  //-->VIDEO5 
double* roulette;
/* + ----------------------------------------------------
   ! VARIABLES UTILES
   +-----------------------------------------------------*/
unsigned chromosome_length; //Tamaño del cromosoma
double crossover_probability; //(--->video 3)
double mutation_probability //(--->video 4)
unsigned max_generations; // video 6
//Video 6 parametros del algortimo
void getParameters() 
{
    printf("\n Numero maximo de generaciones: "); scanf("%d", &max_generations);
    printf("\n Probabilidad de cruza: "); scanf("%lf", &crossover_probability);
    printf("\n Probabilidad de mutaciones: "); scanf("%lf", &mutation_probability);
    
}

/* + ----------------------------------------------------
   ! RESERVAR MEMORIA PARA LAS POBLACIONES
   +-----------------------------------------------------*/
void allocateMemory()       //Hacer una función que va a guardar esas poblaciones
{
    unsigned required_bytes = sizeof(Individual) * POPULATION_SIZE;
    parents = (Individual*) malloc(required_bytes);         // Definir la poblacion de padres
    offspring = (Individual*) malloc(required_bytes);       //Definir la desendencia


    //Vasmos ahora a reservar memoria para los cromosomas de manera itereable
    int i;
    chromosome_length = ceil(log2(( FX_UPPER_BOUND - FX_LOWER_BOUND)*pow(10,PRECISION))); //Definido el largo del cromosoma
    required_bytes = chromosome_length * sizeof(int);
    for(i=0; i<POPULATION_SIZE; i++){
        parents[i].chromosome = (int*) malloc(required_bytes);
        offspring[i].chromosome = (int*) malloc(required_bytes);
        // (VAMOS A RESERVAR MEMORIA PARA LA RULETA --> VIDEO 2)
        roulette = (double*) malloc(sizeof(double)* POPULATIONSIZE);
    }
}
/* + ----------------------------------------------------
   ! FUNCION QUE GENERA LOS ALEATORIOS -> devolver un aleatorio sobre un intervalo
   +-----------------------------------------------------*/
double randomDouble(double a, double b) 
{
    return (b-a)*(rand()/(double)RAND_MAX) + a;
}
/* + ----------------------------------------------------
   ! SIMULADO DE LANZAR UNA MODENA AL AIRE
   +-----------------------------------------------------*/
int flip(double probability)
{
    if(randomDouble(0,1) < probability)
        return 0;
    else
        return 1;
}

/* + ----------------------------------------------------
   ! INICIALIZAR LA PRIMERA GENERACIÓN
   +-----------------------------------------------------*/
void createFirstGeneration()
{
    int i;
    for (i=0; i<POPULATION_SIZE;i++) {
        parents[i].x = RAND_MAX;                                    //Fenotipo
        parents[i].fitness = 0;                                     //Aptitud
        parents[i].parents[0] = parents[i].parents[1] = -1;         //Padres
        parents[i].mutation_place = parents[i].crossover_place = -1;//Mutacion y cruce

        ////Cromosoma//////

        int j ;
        for(j = 0; j<chromosome_length; j++){
            parents[i].chromosome[j] = flip(0.5);
        }
    }
}
/* + ----------------------------------------------------
   ! DECODIFICAR EL GENOTIPO EN BINARIO A REAL Y APLICAR AJUSTE AL RANGO
   +-----------------------------------------------------*/
double binary2real(int* chromosome)
{
    //Converión a  numero binario -> codificación
    int i,
    double aux = 0.0;

    for(i = chromosome_length-1 ; i >=0 ; i--)
        aux += (double) pow(2,chromosome_length - i - 1)    
        //Ajuste al rango
    return FX_LOWER_BOUND + ((aux*(FX_UPPER_BOUND-FX_LOWER_BOUND))/(pow(2,chromosome_length)-1));
}
/* + ----------------------------------------------------
   ! EVALUAR APTITUD DEL INDIVIDUO
   +-----------------------------------------------------*/
void evaluateTargetFunction(Individual* Individual)
{
    individual->x = binary2real(individual->chromosome); // le estamos haciendo la decodificación del ajuste al rango
    individual->fitness = 1/(pow(individual->x, 2)+0.001);       //Asignando su aptitud
}
/* + ----------------------------------------------------
   ! EVALUAR APTITUD DE UNA POBLACIÓN
   +-----------------------------------------------------*/
void evaluatePopulation(Individual* population)
{
    int i;
    for(i = 0; i<POPULATION_SIZE; i++)
        evaluateTargetFunction(&population[i]);
}
/* + ----------------------------------------------------
   ! LLENAR LA RULETA CON LA PROBABILIDAD DE CADA INDIVIDUO PARA SER SELECCIONADO
   +-----------------------------------------------------*/
void updateRoulette(Individual* population)
{
    int i;
    double sum_fitness = 0.0;
    for(i = 0; i<POPULATION_SIZE; i++)
        sum_fitness += population[i].fitness;
    for(i = 0; i<POPULATION_SIZE; i++)
        roulette[i] = poblacion[i].fitness / sum_fitness;
}
/* + ----------------------------------------------------
   ! HACER LA SELECCION POR LA RULETA
   +-----------------------------------------------------*/
unsigned rouletteWheelSelection()
{
    double r = randomDouble(0,1);
    double sum = 0.0;
    int i;
    for (i = 0; sum<r; i++)
        sum += roulette[i];
    return i;
}

/* + ----------------------------------------------------
   ! RECOMBINACION DE LOS PADRES SELECCIONADOS
   +-----------------------------------------------------*/
void crossover(Individual* father, Individual* mother, Individual* child1, Individual* child2)
{
    int i;
    if (flip(crossover_probability)){
        unsigned p = (unsigned) randomDouble(1,chromosome_length-2);
        for(i=p+1; i<chromosome_length,i++){
            child1->chromosome[i] = father->chromosome[i];
            child2->chromosome[i-p-1] = father->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = p;
    } else {
        for(i=0; i<chromosome_length,i++){
            child1->chromosome[i] = father->chromosome[i];
            child1->chromosome[i] = father->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = -1;
    }
}
/* + ----------------------------------------------------
   ! MUTACION
   +-----------------------------------------------------*/
void mutation(Individual* individual)
{
    if(flip(mutation_probability)){
        unsigned p = (unsigned) randomDouble(0, chromosome_length-1);
        individual->chromosome[P] = 1 - individual->chromosome[P];
        individual->mutation_place = p;
    }else{
        individual->mutation_place = -1;
    }
}
/* + ----------------------------------------------------
   ! LEY DEL MAS FUERTE - elitismo--video5)
   +-----------------------------------------------------*/

void elitism(){
    unsigned words_child1 = words_child2 = 0;
    the_best = parents[0];

    for(i=0; i<POPULATION_SIZE;i++){
        if(offspring[i].fitness < offspring[words_child1].fitness)
            words_child1 = i;
        else if(offspring[i].fitness < offspring[words_child2].fitness)
            words_child2 = i;

        if(parents[i].fitness > the_best.fitness.fitness)
            the_best = parents[i];
    }

    offspring[words_child1] = the_best;
    offspring[words_child2] = the_best;
    
}



/* + ----------------------------------------------------
   ! FUNCION PRINCIPAL
   +-----------------------------------------------------*/
   /* + ----------------------------------------------------
   ! FUNCION PRINCIPAL
   +-----------------------------------------------------*/
int main ()
{
    /* + ----------------------------------------------------
   ! PREPARACION
   +-----------------------------------------------------*/
    getParameters()                     //video6
    srand(long)time(NULL);              //Generación de aleatorios

    allocateMemory();
    createFirstGeneration();
    evaluatePopulation(parents);
    updateRoulette(parents);            //llenar la ruleta video 6
    int generation;

    for (generation = 0; generation<max_generations; generation++){
        /* + ----------------------------------------------------
        ! PROCESO DE GENERACIÓN video 6
        +-----------------------------------------------------*/
        for(i = 0; i < POPULATION_SIZE-1;i+=2){ 
            unsigned selected_father = rouletteWheelSelection();     //seleccionar individuos  video 6
            unsigned selected_mother = rouletteWheelSelection();     //seleccionar individuos  video 6
            crossover(&parents[selected_father], &parents[selected_mother], &offspring[i], &offspring[i+1]);   //Los cruzamos video 6
            mutation(&offspring[i]);                                //Mutamos a los hijos video 6
            mutation(&offspring[i+1]);                              //Mutamos a los hijos video 6 
        }
        evaluatePopulation(offspring);   // Ya tenemos el valor de aptitud de los hijos, el paso siguiente es aplicar el elitismo
        elitism();                       //Aplicamos elitimos
        parents = offspring;            // Ya los hijos sustituyen a los padres video 6  solo para que funcione
        updateRoulette(parents);        // Actualizamos los valores de la ruleta video 6 solo para que funcione
    }   
    return 0;
}



//Notas : El video 6 es el correspondiente a la clase 28 de Udemy, explicación de todo minuto 5:37