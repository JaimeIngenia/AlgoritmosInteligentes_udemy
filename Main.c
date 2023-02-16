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
double mutation_probability; //(--->video 4)
unsigned max_generations; // video 6
unsigned selected_father;  ///CORRECION VIDEO 7 
unsigned selected_mother;  ///CORRECION VIDEO 7


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
    //required_bytes = chromosome_length * sizeof(int); //CORRECION VIDEO 7

    for(i=0; i<POPULATION_SIZE; i++){
        parents[i].chromosome = (int*) calloc(chromosome_length, sizeof(int));
        offspring[i].chromosome = (int*) calloc(chromosome_length, sizeof(int));
        parents[i].x = offspring[i].x = RAND_MAX;  //correcion video 7
        parents[i].fitness = offspring[i].fitness = 0;         //correcion video 7
        parents[i].parents[0] = offspring[i].parents[0] = -1;         //correcion video 7
        parents[i].parents[1] = offspring[i].parents[1] = -1;         //correcion video 7
        parents[i].crossover_place = offspring[i].crossover_place = -1;         //correcion video 7
        parents[i].mutation_place = offspring[i].mutation_place = -1;         //correcion video 7
    }
    the_best.chromosome = (int*) calloc(chromosome_length, sizeof(int));   //correcion video 7
    // (VAMOS A RESERVAR MEMORIA PARA LA RULETA --> VIDEO 2) correcion video 7
    roulette = (double*) malloc(sizeof(double)* POPULATION_SIZE);
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
    if(randomDouble(0,1) <= probability)
        return 1;  //correcion video 7
    else
        return 0;  //correcion video 7
}

/* + ----------------------------------------------------
   ! INICIALIZAR LA PRIMERA GENERACIÓN
   +-----------------------------------------------------*/
void createFirstGeneration()
{
    int i, j ;
    for (i=0; i<POPULATION_SIZE;i++)          //correcion video 7
        for(j = 0; j<chromosome_length; j++)     //correcion video 7
            parents[i].chromosome[j] = flip(0.5); //correcion video 7  // generacion del cromosoma aleatorio
}
/* + ----------------------------------------------------
   ! DECODIFICAR EL GENOTIPO EN BINARIO A REAL Y APLICAR AJUSTE AL RANGO
   +-----------------------------------------------------*/
double binary2real(int* chromosome)
{
    //Converión a  numero binario -> codificación
    int i;
    double aux = 0.0;
    for(i = chromosome_length-1 ; i >=0 ; i--){
        if(chromosome[i] == 1) //correcion video 7                                                                                                                                                                                                                                                                                                                                                                                                      
            aux += (double) pow(2,chromosome_length - i - 1); 
        //Ajuste al rango
    }
    return FX_LOWER_BOUND + ((aux*(FX_UPPER_BOUND-FX_LOWER_BOUND))/(pow(2,chromosome_length)-1));
}
/* + ----------------------------------------------------
   ! EVALUAR APTITUD DEL INDIVIDUO
   +-----------------------------------------------------*/
void evaluateTargetFunction(Individual* individual)
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
        roulette[i] = population[i].fitness / sum_fitness;
}
/* + ----------------------------------------------------
   ! HACER LA SELECCION POR LA RULETA
   +-----------------------------------------------------*/
unsigned rouletteWheelSelection()
{
    double r = randomDouble(0,1);
    double sum = 0.0;
    int i, current_individual; //CORRECCIO VIDEO 7
    for (i = POPULATION_SIZE; sum<r; i++){
        current_individual = i % POPULATION_SIZE; //CORRECCIO VIDEO 7
        sum += roulette[current_individual];
    }
    return current_individual;    //CORRECCIO VIDEO 7
}

/* + ----------------------------------------------------
   ! RECOMBINACION DE LOS PADRES SELECCIONADOS
   +-----------------------------------------------------*/
void crossover(Individual* father, Individual* mother, Individual* child1, Individual* child2)
{
    int i;
    if (flip(crossover_probability)){
        unsigned p = (unsigned) randomDouble(1,chromosome_length-2);
        for(i=0; i<p;i++){
            child1->chromosome[i] = father->chromosome[i];
            child2->chromosome[p+i] = mother->chromosome[i];
        }
        for(i=p+1; i<chromosome_length;i++){
            child1->chromosome[i] = mother->chromosome[i];
            child2->chromosome[i-p-1] = father->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = p;
        child1->parents[0] = child2->parents[0] = selected_father +1;   //COORECOPN VIDEO 7
        child1->parents[1] = child2->parents[1] = selected_mother +1;   //COORECOPN VIDEO 7
    } else {
        for(i=0; i<chromosome_length;i++){
            child1->chromosome[i] = father->chromosome[i];
            child2->chromosome[i] = mother->chromosome[i];
        }
        child1->crossover_place = child2->crossover_place = -1;
        child1->parents[0] = child2->parents[0] = 0;   //COORECOPN VIDEO 7
        child1->parents[1] = child2->parents[1] = 0;   //COORECOPN VIDEO 7
    }
}
/* + ----------------------------------------------------
   ! MUTACION
   +-----------------------------------------------------*/
void mutation(Individual* individual)
{
    if(flip(mutation_probability)){
        unsigned p = (unsigned) randomDouble(0, chromosome_length-1);
        individual->chromosome[p] = 1 - individual->chromosome[p];
        individual->mutation_place = p;
    }else{
        individual->mutation_place = -1;
    }
}
/* + ----------------------------------------------------
   ! LEY DEL MAS FUERTE - elitismo--video5)
   +-----------------------------------------------------*/

void elitism(){
    unsigned best_parent;
    unsigned words_child1, words_child2 ;
    int i;

    best_parent = words_child1 = words_child2 = 0;
    

    for(i=0; i<POPULATION_SIZE;i++){
        if(offspring[i].fitness < offspring[words_child1].fitness)
            words_child1 = i;
        else if(offspring[i].fitness < offspring[words_child2].fitness)
            words_child2 = i;

        if(parents[i].fitness > parents[best_parent].fitness)
            best_parent = i;
    }

    offspring[words_child1] = parents[best_parent];
    offspring[words_child2] = parents[best_parent];
    
}
/* + ----------------------------------------------------
   ! IMPRIMIR UN CROMOSOMA CON INDICADORES   video 7
   +-----------------------------------------------------*/
void printChromosome(Individual* individual)
{
    int i;
    for(i=0; i<chromosome_length; i++)
        if(i == individual->mutation_place) printf("(");
        printf("%d ", individual->chromosome[i]);
        if(i == individual->mutation_place) printf(")");
        if(i == individual->crossover_place) printf("/");
}

/* + ----------------------------------------------------
   ! MOSTRAR INFORMACION DE UNA POBLACIÓN
   +-----------------------------------------------------*/
void printPopulacionDetail(Individual* population)
{
    int i;
    printf("\n\n-------------------------------\n");
    printf(" #\tChromosome\tx\tFitness\tParents");
    printf("\n\n-------------------------------\n");
    for(i=0; i<POPULATION_SIZE; i++){
        printf("\n%03d ", i+1);
        printChromosome(&population[i]);
        printf(" %.3f\t%.3f\t(%d,%d)", population[i].x, population[i].fitness, 
                                       population[i].parents[0],population[i].parents[1]);
    }
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
    getParameters();                    //video6
    srand((long)time(NULL));              //Generación de aleatorios

    allocateMemory();
    createFirstGeneration();
    evaluatePopulation(parents);
    Individual* temp_helper;          //Correción video 7
    //updateRoulette(parents);            //llenar la ruleta video 6
    int generation,i;

    for (generation = 0; generation<max_generations; generation++){
        updateRoulette(parents); 
        printPopulacionDetail(parents);
        /* + ----------------------------------------------------
        ! PROCESO DE GENERACIÓN video 6
        +-----------------------------------------------------*/
        for(i = 0; i < POPULATION_SIZE-1;i+=2){ 
            selected_father = rouletteWheelSelection();     //seleccionar individuos  video 6
            selected_mother = rouletteWheelSelection();     //seleccionar individuos  video 6
            crossover(&parents[selected_father], &parents[selected_mother], &offspring[i], &offspring[i+1]);   //Los cruzamos video 6
            mutation(&offspring[i]);                                //Mutamos a los hijos video 6
            mutation(&offspring[i+1]);                              //Mutamos a los hijos video 6 
            evaluateTargetFunction(&offspring[i]);    // correcion video 7
            evaluateTargetFunction(&offspring[i+1]);    // correcion video 7
        }
        //evaluatePopulation(offspring);   // Ya tenemos el valor de aptitud de los hijos, el paso siguiente es aplicar el elitismo
        elitism();                       //Aplicamos elitimos
        temp_helper = parents;
        parents = offspring;            // Ya los hijos sustituyen a los padres video 6  solo para que funcione
        offspring = temp_helper;
        //updateRoulette(parents);        // Actualizamos los valores de la ruleta video 6 solo para que funcione
        printf("\n\n\tTermino con exito generación %d\n", generation+1);
    }   
    free(parents);
    free(offspring);
    free(roulette);
    return 0;
}



//Notas : El video 6 es el correspondiente a la clase 28 de Udemy, explicación de todo minuto 5:37