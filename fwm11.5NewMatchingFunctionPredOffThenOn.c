/* Blakes Food Web model
 copyright Blake Stevison, 2014, all rights reserved
 version 1.6, FoodWebModel.c */


#include <stdio.h>      // Packages to Install For Model //
#include <time.h>       // Packages to Install For Model //
#include <math.h>       // Packages to Install For Model //
#include <stdlib.h>     // Packages to Install For Model //
#include <string.h>     // Packages to Install For Model //
#include <curses.h>     // Packages to Install For Model //
#include <gsl_rng.h>    // gnu scientific library //
#include <gsl_randist.h>// gnu scientific library //
#include <getopt.h>  // SMF added 10.14.14
#include <gsl_statistics_double.h> // package to run correlation

//////////////////////////////////////////////////////////////////////////////////////

// parameter defaults initialization //

#define deltaT_DEFAULT 0.1 // deltaT is the time scaler to change when things update

#define NS 500 // default number of species in whole system

#define T_TOTAL 5000 // default total number of time steps

#define K_DEFAULT 5000.0 // primary producer carrying capacity, was 5000

#define rPP_DEFAULT 1.05 // was 1.2 // primary producer intrisic rate of increase

//////////////////////////////////////////////////////////////////////////////////////

#define PROPPRIMS_DEFAULT 0.6 // proportion of primary producers in the food web

#define PROPHERBS_DEFAULT 0.4 // proportion of herbivores in the food web

#define PROPPREDS_DEFAULT 0.0 // proportion of predators in the food web

//////////////////////////////////////////////////////////////////////////////////////

#define MEANINITConsPOPSIZE_DEFAULT 2500.0 // mean initial population sizes

#define SDINITConsPOPSIZE_DEFAULT 800.0 // standard deviation initial pop sizes

#define MEANINITPrimPOPSIZE_DEFAULT 2500.0 // mean initial population sizes

#define SDINITPrimPOPSIZE_DEFAULT 800.0 // std dev of initial population sizes

//////////////////////////////////////////////////////////////////////////////////////

#define AttackRate_DEFAULT 0.2 // proportion of time when start collecting GainLossAll

#define BASALMR_DEFAULT 0.05 // Basal metabolic rate for consumers

#define EFFORTINCR_DEFAULT 0.01 // foraging effort incriment looped over

#define THETA_DEFAULT 20.0 //54.598 // Theta = diet breadth coefficient for consumers

//////////////////////////////////////////////////////////////////////////////////////

#define ADDEDSPECIES_DEFAULT 0 // number of species added at certain time steps

#define REMOVEDSPECIES_DEFAULT 0 // number of species removed at cetain time steps

#define GAINSSCALER_DEFAULT 0.1 // divide gains by to scale populaiton dynamics

#define LOWSPACE_DEFAULT 0.0 // used to assign producers/consumers in nutrient space

#define HIGHSPACE_DEFAULT 0.45 // assign producers and consumers in nutrient space

#define PARAMSTEPS_DEFAULT 1 // number of loops over the parameter being changed

#define PRINTTIME_DEFAULT 10 // print to file every PRINTTIME time steps

#define StartCollectTime_DEFAULT 0.6 // proportion of time when collecting GainLossAll

#define SigCorr_DEFAULT 0.159 // significant correlation cutt of to define Sig

#define ChangeCritical_DEFAULT 1.0 // critical chagne amount to stop parameter run

//////////////////////////////////////////////////////////////////////////////////////

int KillProducerTime = 25000;

// parameters that are global variables
double deltaT = deltaT_DEFAULT ;
//int T_TOTAL = T_TOTAL_DEFAULT;
int Kcap = K_DEFAULT;
double rPP = rPP_DEFAULT;
double PROPPRIMS = PROPPRIMS_DEFAULT;
double PROPHERBS = PROPHERBS_DEFAULT;
double PROPPREDS = PROPPREDS_DEFAULT;
int MEANINITPrimPOPSIZE = MEANINITPrimPOPSIZE_DEFAULT;
int SDINITPrimPOPSIZE = SDINITPrimPOPSIZE_DEFAULT;
int MEANINITConsPOPSIZE = MEANINITConsPOPSIZE_DEFAULT;
int SDINITConsPOPSIZE = SDINITConsPOPSIZE_DEFAULT;
double BASALMR = BASALMR_DEFAULT;
double EFFORTINCR = EFFORTINCR_DEFAULT;
double THETA = THETA_DEFAULT;
int ADDEDSPECIES = ADDEDSPECIES_DEFAULT;
int REMOVEDSPECIES = REMOVEDSPECIES_DEFAULT;
double GAINSSCALER = GAINSSCALER_DEFAULT;
double LOWSPACE = LOWSPACE_DEFAULT;
double HIGHSPACE = HIGHSPACE_DEFAULT;
int PARAMSTEPS = PARAMSTEPS_DEFAULT;
int PRINTTIME = PRINTTIME_DEFAULT;
double StartCollectTime = StartCollectTime_DEFAULT;
double AttackRate = AttackRate_DEFAULT;
double SigCorr = SigCorr_DEFAULT;
double ChangeCritical = ChangeCritical_DEFAULT;

int MinTimeStable = (T_TOTAL-1);


int TrophicLevel[NS];
// globals
int prim[NS]; // True or False if primary producer, prim=1=primary prod, prim=0=consumer
int pred[NS];
double percentC[NS]; // percent of organism's  cabron; used to be x1 in mathematica code
double percentN[NS]; // percent of orgainsm's nitrogen; used to be x2 in mathematica code
double percentP[NS]; // percent of orgainsm's Phosphorus; used to be x2 in mathematica
double Cent[NS];
double PopSize[NS]; // population size of all species
double PopSizePrevious[NS]; // place holder for the previeous population size
double K[NS];
int aPrey[NS]; //the a prey item that each species eats
int bPrey[NS]; //the b prey item that each species eats
int cPrey[NS]; //the c prey item that each species eats
double Edist[NS];
double EdistOpt[NS];
double OptEffa[NS];
double OptEffb[NS];
double OptEffc[NS];
double TotalKilled[NS];
double Gains[NS];
int t;
int chainLength[NS];
int numCons = 0;
int numPrims = 0;
int numPrimsForced;
int numConsForced;
int numPredsForced;

// Initializing percent element vectors
double percentC_aPrey[NS];
double percentC_bPrey[NS];
double percentC_cPrey[NS];
double percentN_aPrey[NS];
double percentN_bPrey[NS];
double percentN_cPrey[NS];
double percentP_aPrey[NS];
double percentP_bPrey[NS];
double percentP_cPrey[NS];

// Calculating percent C,N,P of diet from prey items and efforts, x1E in mathematica
double percentC_Diet[NS];
double percentN_Diet[NS];
double percentP_Diet[NS];

// normalizing vectors of diet compositoins
double percentC_Diet_Norm[NS];
double percentN_Diet_Norm[NS];
double percentP_Diet_Norm[NS];

// from OptForgFunction
double Edist[NS];
double EdistOpt[NS];
double rawgains[NS];
double Optrawgains[NS];

// GainLossFunction Variables
double GainLoss[NS][NS];
double GainLossAll[NS][NS];
double EffortMatrix[NS][NS];

// from MainLoop
double MeanPopSize;
double MeanConsPopSize;
double MeanPrimPopSize;
int SpeciesAlive[NS];
int Extinctions[T_TOTAL];
int numConsAlivePrev;
int numPrimsAlivePrev;
// Triangles
int Triangle[NS];
double MeanNonTrisPopSize;
double MeanTrueTrisPopSize;

// tracking stabalization
double PopSizeChange[NS];
double MeanPopSizeChange;
double MeanPopSizeChange;

double TotalPopSizeChange;

// from initialization and choose functions
int numPrimsForced,numHerbsForced,numPredsForced;

int Stable_Trigger;

// globals for GSL random number generation
const gsl_rng_type *T;   /*generator type*/
gsl_rng *r;             /*rng instance*/


#define randu() gsl_rng_uniform(r) // this calls for uniform random number

// function declarations
double EucDistFunction(double x1, double y1, double z1, double x2, double y2, double z2);

void InitializationFunction(unsigned long int SeedSet);

void cleanUpAndFinishFunction(void);

int randInt(int min, int max);

void chooseFunction(int nChoices);

void OptForagingFunction(void);

void GainLossFunction(void);

void PopulationUpdateFunction(void);

void usage(char *s); // SMF added 10.14.14

void calcChainLength(void);

void PrintParams(void);

void DropExtremeSpecies(void);

void UnChangingInitializationFunction(void);

// declare files
FILE * TimeFile;
FILE * InitFile;
FILE * OptForgFile;
FILE * GainLossFile;
FILE * GainLossAllFile;
FILE * PopUpdateFile;
FILE * ConsInitFile;
FILE * PrimInitFile;
FILE * ConsPopUpdateFile;
FILE * PrimPopUpdateFile;
FILE * ConsOptForgFile;
FILE * PrimOptForgFile;
FILE * TrophicFile;
FILE * EffortMatrixFile;
FILE * ParamFile;
FILE * ParamPopSizeFile;
FILE * TimeFileAllPops;
FILE * ppTimeFileAllPops;
FILE * consTimeFileAllPops;
FILE * ExtinctionLogFile;
FILE * SpeciesAliveParamFile;
FILE * LastTimeFile;
FILE * TriangleFileInit;
FILE * TriangleFile;
FILE * TriTimeFile;


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// main loop and beginnning of action
int main(int argc, char *argv[]) {
    
    
    TimeFile = fopen("TimeFile.csv", "w");
    InitFile = fopen("InitFile.csv", "w");
    OptForgFile = fopen("OptForgFile.csv", "w");
    GainLossFile = fopen("GainLossFile.csv","w");
    GainLossAllFile = fopen("GainLossAllFile.csv","w");
    PopUpdateFile = fopen("PopUpdateFile.csv","w");
    ConsInitFile = fopen("ConsInitFile.csv","w");
    PrimInitFile = fopen("PrimInitFile.csv","w");
    ConsPopUpdateFile = fopen("ConsPopUpdateFile.csv","w");
    PrimPopUpdateFile = fopen("PrimPopUpdateFile.csv","w");
    ConsOptForgFile = fopen("ConsOptForgFile.csv","w");
    PrimOptForgFile = fopen("PrimOptForgFile.csv","w");
    TrophicFile = fopen("TrophicFile.csv","w");
    EffortMatrixFile = fopen("EffortMatrixFile.csv","w");
    ParamFile = fopen("ParamFile.csv","w");
    ParamPopSizeFile = fopen("ParamPopSizeFile.csv","w");
    TimeFileAllPops = fopen("TimeFileAllPops.csv","w");
    ppTimeFileAllPops = fopen("ppTimeFileAllPops.csv","w");
    consTimeFileAllPops = fopen("consTimeFileAllPops.csv","w");
    ExtinctionLogFile = fopen("ExtinctionLogFile.csv","w");
    SpeciesAliveParamFile = fopen("SpeciesAliveParamFile.csv","w");
    LastTimeFile = fopen("LastTimeFile.csv","w");
    TriangleFile = fopen("TriangleFile.csv","w");
    TriTimeFile = fopen("TriTimeFile.csv","w");
    
    int i;
    int ch; // SMF added 10.14.14
    char *progname = argv[0]; // SMF added 10.14.14
    double TotalPopSize;
    double TotalConsPopSize;
    double TotalPrimPopSize, NSdubs;
    
    //TriangleMeanCalcs
    double TotalTrueTrisPopSize;
    double TotalNonTrisPopSize;
    
    // correlation data
    double corr;
    double corrPrims;
    double corrCons;
    int Sig;
    int SigPrims;
    int SigCons;
    double var;
    double varPrims;
    double varCons;
    unsigned long int SeedSet = 2;
    
    /* beginning of a block of new stuff added by SMF 10.14.14 */
    while ((ch = getopt(argc, argv, "d:K:P:H:D:R:t:A:M:S:?")) != -1) {
        // program should now accept -D -K and -T "flags" on command line
        // ":" after each letter means that the flag needs an argument, i.e., a number in the appropriate format
        switch (ch) {
                case 'd':
                deltaT = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'K':
                Kcap = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'P':
                PROPPRIMS = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'H':
                PROPHERBS = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'D':
                PROPPREDS = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'R':
                rPP = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 't':
                THETA = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'A':
                AttackRate = strtod(optarg, (char **)NULL); // read in as double
                break;
                case 'M':
                BASALMR = strtod(optarg, (char **)NULL); // read in as double
                break;
            case 'S':
                SeedSet = strtod(optarg, (char **)NULL); // read in as double
                break;
                case '?':
            default:
                usage(progname); // a function for human readers to be able to see what each thing does
                exit(-1);
        }
    }
    /* end of a block of new stuff added by SMF 10.14.14 */
    
    // sam suggests putting RNG initialization here and removing it from
    // the "InitializationFunction()"
    
    
    //for (Community = 0; Community < 1; Community++) {
    
    // initialize random number generator
    gsl_rng_env_setup(); // set up the environment variables for the RNG
    T = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    
    // Function Calls
    UnChangingInitializationFunction();
    
    // initialize all compositions and populaiton size
    InitializationFunction(SeedSet);
    
    // choose prey items
    chooseFunction(3);
    
    NSdubs = (double) NS;
    
    t=0;
    
    Stable_Trigger = 1;
    
    int numConsAlive;
    int numPrimsAlive;
    int numNonTris;
    int numTrueTris;
    
    
    while (t<T_TOTAL && Stable_Trigger == 1) {
        
        // reset for each time step
        TotalPopSize = 0.0;
        TotalPrimPopSize = 0.0;
        TotalConsPopSize = 0.0;
        TotalPopSizeChange = 0.0;
        TotalTrueTrisPopSize = 0.0;
        TotalNonTrisPopSize = 0.0;
        
        // find optimal foraging efforts
        OptForagingFunction();
        
        // find raw gains and losses
        GainLossFunction();
        
        // calcualte population dynamics, given foraging efforts
        PopulationUpdateFunction();
        
        // SMF QUESTION: Why is all this stuff here instead of in the population update function?
        
        // calculating mena population size each time step
        if (t % PRINTTIME == 0.0) {
            
            fprintf(TimeFileAllPops,"%i",t);
            fprintf(ppTimeFileAllPops,"%i",t);
            fprintf(consTimeFileAllPops,"%i",t);
            
        }
        
        numConsAlivePrev = numConsAlive;
        numPrimsAlivePrev = numPrimsAlive;
        
        numConsAlive = 0;
        numPrimsAlive = 0;
        numNonTris = 0;
        numTrueTris = 0;
        
        for (i=0; i<NS; i++) {
            
            if(Triangle[i]) {
            
            printf("%i\n", Triangle[i]);
            
            }
            
            if (PopSize[i] > 2.0) {
                
                if ( prim[i])
                numPrimsAlive++;
                else
                numConsAlive++;
            }
                
                if ( Triangle[i])
                    numTrueTris++;
            
                if ( !Triangle[i] && !prim[i])
                    numNonTris++;
            
            TotalPopSize = TotalPopSize + PopSize[i];
            
            if (t % PRINTTIME == 0.0) {
                
                fprintf(TimeFileAllPops,",%f",PopSize[i]);
                
            }
            
            if (prim[i] == 0) {
                
                TotalConsPopSize = TotalConsPopSize + PopSize[i];
                
                if (t % PRINTTIME == 0.0) {
                    
                    fprintf(consTimeFileAllPops,",%f",PopSize[i]);
                    
                }
            }
            if (prim[i] == 1) {
                
                TotalPrimPopSize = TotalPrimPopSize + PopSize[i];
                
                if (t % PRINTTIME == 0.0) {
                    
                    fprintf(ppTimeFileAllPops,",%f",PopSize[i]);
                    
                }
            }
            
            // Triangle Mean PopSize Calculations
            //TotalPopSize = TotalPopSize + PopSize[i];
            
            if (Triangle[i] == 0 && !prim[i]) {
                
                TotalNonTrisPopSize = TotalNonTrisPopSize + PopSize[i];
                
                 //printf("\n Got to NonTris \n");
            }
            
            if (Triangle[i] == 1 && !prim[i]) {
                
                TotalTrueTrisPopSize = TotalTrueTrisPopSize + PopSize[i];
                
                printf("\n Got to TrueTris \n");
            }
            
        } // End of NS loop
        
        if (t % PRINTTIME == 0.0) {
            
            fprintf(TimeFileAllPops,"\n");
            fprintf(ppTimeFileAllPops,"\n");
            fprintf(consTimeFileAllPops,"\n");
            
        }
        
        // track stability of the system
        /*MeanPopSizeChange = TotalPopSizeChange / (numPrimsAlive + numConsAlive);
         
         if (MeanPopSizeChange < ChangeCritical) {
         
         t = T_TOTAL;
         }*/
        
        
        MeanPopSize =  TotalPopSize / (((double) numConsAlive) + ((double) numPrimsAlive));
        MeanConsPopSize = TotalConsPopSize / ((double) numConsAlive);
        MeanPrimPopSize = TotalPrimPopSize / ((double) numPrimsAlive);
        // Triangle Means
        MeanTrueTrisPopSize = TotalTrueTrisPopSize / ((double) numTrueTris);
        MeanNonTrisPopSize = TotalNonTrisPopSize / ((double) numNonTris);
        
        // print that shit to a fucking file
        
        if (t % PRINTTIME == 0.0) {
            
            fprintf(TimeFile, "%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i\n", t, MeanPopSize,MeanPrimPopSize,MeanConsPopSize,TotalPopSize,TotalPrimPopSize,TotalConsPopSize,PROPPRIMS,PROPHERBS,PROPPREDS,MeanTrueTrisPopSize,MeanNonTrisPopSize,numNonTris,numTrueTris);
            
        }
        
        //printf("%i %f %f %f \n", t, MeanPopSize,MeanPrimPopSize,MeanConsPopSize);
        
        
        if (t == KillProducerTime) {
            
            //DropExtremeSpecies();
        }
        
        Extinctions[t] = (numPrimsAlive + numConsAlive) - (numPrimsAlivePrev + numConsAlivePrev);
        
        /*if (t > 100000 && (Extinctions[t] - Extinctions[t-1] - Extinctions[t-2] - Extinctions[t-3] - Extinctions[t-4] - Extinctions[t-5] - Extinctions[t-6] - Extinctions[t-7] - Extinctions[t-8] - Extinctions[t-9] - Extinctions[t-10]) == 0) {
         
         t = T_TOTAL;
         }*/
        
        t++;
    }
    
    // find correlation of cent and popsize
    //double CentAlive[NS];
    //double ppCentAlive[numPrimsAlive];
    //double consCentAlive[numConsAlive];
    
    double CentAlive[NS];
    double ppCentAlive[numPrimsAlive];
    double consCentAlive[numConsAlive];
    int TotalPrimsAlive;
    int TotalConsAlive;
    
    int countPrim = 0;
    int countCons = 0;
    
    for (i=0; i<NS; i++) {
        
        if (PopSize[i] > 2.0) {
            
            CentAlive[i] = Cent[i];
            
            if (prim[i]) {
                
                ppCentAlive[countPrim] = Cent[i];
                
                countPrim++;
                
            }
            
            else {
                
                consCentAlive[countCons] = Cent[i];
                
                countCons++;
            }
        }
    }
    
    TotalPrimsAlive = countPrim;
    TotalConsAlive =  countCons;
    
    corr = gsl_stats_correlation (&Cent[0], 1, &PopSize[0], 1, NS);
    
    var = gsl_stats_variance (&Cent[0], 1, NS);
    
    //corrPrims = gsl_stats_correlation (&ppCentAlive[0], 1, &PopSize[0], 1, NS);
    
    varPrims = gsl_stats_variance (&ppCentAlive[0], 1, numPrimsAlive);
    
    //corrCons = gsl_stats_correlation (&Cent[0], 1, &PopSize[0], 1, NS);
    
    varCons = gsl_stats_variance (&consCentAlive[0], 1, numConsAlive);
    
    // For the whole system //
    if (fabs(corr) > SigCorr) {
        Sig = 1;
    }
    else Sig=0;
    
    // correlations for primary producers
    corrPrims = gsl_stats_correlation (&Cent[0], 1, &PopSize[0], 1, numPrimsForced);
    
    //varPrims = gsl_stats_variance (&Cent[0], 1, numPrimsForced);
    
    if (fabs(corrPrims) > SigCorr) {
        SigPrims = 1;
    }
    else SigPrims=0;
    
    // correlation for consumers
    corrCons = gsl_stats_correlation (&Cent[numPrimsForced], 1, &PopSize[numPrimsForced], 1, (NS-numPrimsForced));
    
    //varCons = gsl_stats_variance (&Cent[numPrimsForced], 1, (NS - numPrimsForced));
    
    if (fabs(corrCons) > SigCorr) {
        SigCons = 1;
    }
    else SigCons=0;
    
    for (i=0; i<NS; i++) {
        
        if (PopSize[i] > 2.0) {
            
            SpeciesAlive[i] = 1;
        }
        else SpeciesAlive[i] = 0;
        
        fprintf(SpeciesAliveParamFile,"%i,%f,%f,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%f,%f,%f,%i,%i\n",i,Edist[i],Cent[i],SpeciesAlive[i],prim[i],pred[i],percentC[i],percentN[i],percentP[i],PopSize[i],corr,corrCons,corrPrims,var,varPrims,varCons,PROPPRIMS,PROPHERBS,PROPPREDS,THETA,Kcap,rPP,AttackRate,BASALMR,T_TOTAL,NS);
    }
    
    for ( i=0; i<NS; i++ ) {
        
        float p0x[NS];
        float p1x[NS];
        float p2x[NS];
        
        float p0y[NS];
        float p1y[NS];
        float p2y[NS];
        
        float Area[NS];
        float s[NS];
        float tt[NS];
        float p[NS];
        
        float px[NS];
        float py[NS];
        
        px[i] = percentC[i];
        py[i] = percentN[i];
        
        // percent carbon composition of prey
        percentC_aPrey[i] = percentC[aPrey[i]];
        percentC_bPrey[i] = percentC[bPrey[i]];
        percentC_cPrey[i] = percentC[cPrey[i]];
        // for traingel algorithm
        p0x[i] = percentC[aPrey[i]];
        p1x[i] = percentC[bPrey[i]];
        p2x[i] = percentC[cPrey[i]];
        
        // percent nitrogen composition of prey
        percentN_aPrey[i] = percentN[aPrey[i]];
        percentN_bPrey[i] = percentN[bPrey[i]];
        percentN_cPrey[i] = percentN[cPrey[i]];
        // for triangle algorithm
        p0y[i] = percentN[aPrey[i]];
        p1y[i] = percentN[bPrey[i]];
        p2y[i] = percentN[cPrey[i]];
        
        // percent phosphorus composition of prey
        percentP_aPrey[i] = percentP[aPrey[i]];
        percentP_bPrey[i] = percentP[bPrey[i]];
        percentP_cPrey[i] = percentP[cPrey[i]];
        
        // determines if the composition is in or outside the triangle
        Area[i] = 0.5 * ( (-p1y[i] * p2x[i]) + (p0y[i] * (-p1x[i] + p2x[i] )) + (p0x[i] * (p1y[i]  - p2y[i] )) + (p1x[i] * p2y[i]) );
        //printf("\narea = %f, %f\n", Area, (p0y * (-p1x  + p2x )));
        
        s[i]  = (1.0 / (2.0 * Area[i])) * (p0y[i] * p2x[i]  - p0x[i] * p2y[i] + (p2y[i]  - p0y[i] ) * px[i]  + (p0x[i]  - p2x[i] ) * py[i] );
        tt[i]  = (1.0 / (2.0 * Area[i])) * (p0x[i] * p1y[i]  - p0y[i] * p1x[i] + (p0y[i]  - p1y[i] ) * px[i]  + (p1x[i]  - p0x[i] ) * py[i] );
        
        if (s[i]  >= 0 && s[i]  <= 1 && tt[i] >= 0 && tt[i]  <= 1 && (s[i]  + tt[i] ) <= 1) {
            
            Triangle[i]  = 1;
        }
        
        else {
            
            Triangle[i]  = 0;
        }
        
        fprintf(TriangleFile,"%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i,%f,%f\n", i, Triangle[i],s[i],tt[i],p0x[i],p1x[i],p2x[i],p0y[i],p1y[i],p2y[i],aPrey[i],bPrey[i],cPrey[i],Cent[i],PopSize[i]);
    }
    
    // Calculate Mean Tri Pops
    numNonTris = 0;
    numTrueTris = 0;
    
    for (i=0; i<NS; i++) {
        
        
        if ( Triangle[i])
            numTrueTris++;
        
        if ( !Triangle[i] && !prim[i])
            numNonTris++;

        // Triangle Mean PopSize Calculations
        
        if (Triangle[i] == 0 && !prim[i]) {
            
            TotalNonTrisPopSize = TotalNonTrisPopSize + PopSize[i];
            
            //printf("\n Got to NonTris \n");
        }
        
        if (Triangle[i] == 1 && !prim[i]) {
            
            TotalTrueTrisPopSize = TotalTrueTrisPopSize + PopSize[i];
            
            printf("\n Got to TrueTris \n");
        }
        
    } // End of NS loop
    
    // Triangle Means
    MeanTrueTrisPopSize = TotalTrueTrisPopSize / ((double) numTrueTris);
    MeanNonTrisPopSize = TotalNonTrisPopSize / ((double) numNonTris);
    
    // print that shit to a fucking file
    
    if (t % PRINTTIME == 0.0) {
        
        fprintf(TriTimeFile, "%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i\n", t, MeanPopSize,MeanPrimPopSize,MeanConsPopSize,TotalPopSize,TotalPrimPopSize,TotalConsPopSize,PROPPRIMS,PROPHERBS,PROPPREDS,MeanTrueTrisPopSize,MeanNonTrisPopSize,numNonTris,numTrueTris);
        
    }
    
    //}
    
    //function that calculates trophic web
    calcChainLength();
    
    // print parameter values to file
    PrintParams();
    
    // Close all the files
    fclose(TimeFile);
    fclose(TriTimeFile);
    fclose(InitFile);
    fclose(OptForgFile);
    fclose(GainLossFile);
    fclose(GainLossAllFile);
    fclose(PopUpdateFile);
    fclose(ConsInitFile);
    fclose(PrimInitFile);
    fclose(ConsPopUpdateFile);
    fclose(PrimPopUpdateFile);
    fclose(ConsOptForgFile);
    fclose(PrimOptForgFile);
    fclose(TrophicFile);
    fclose(EffortMatrixFile);
    fclose(ParamFile);
    fclose(ParamPopSizeFile);
    fclose(TimeFileAllPops);
    fclose(ppTimeFileAllPops);
    fclose(consTimeFileAllPops);
    fclose(ExtinctionLogFile);
    fclose(SpeciesAliveParamFile);
    fclose(LastTimeFile);
    fclose(TriangleFile);
    
    return 0;
    
    // all done, clean up
    //cleanUpAndFinishFunction();
    
} // end of main loop

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// defining all functions bellow // // defining all functions bellow // // defining all functions bellow //

void cleanUpAndFinishFunction(void)

{
    // use this function for any final prints and checks you want to do at the end of a simulation
    gsl_rng_free(r); // this frees the memory that was dynamically allocated to the RNG
}

//////////////////////////////////////////////////////////////////////////////////////////

void UnChangingInitializationFunction(void)
{
    
    // initialize random number generator
    gsl_rng_env_setup(); // set up the environment variables for the RNG
    T = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    
    int i;
    double dumC, dumN, dumP, sumdum;
    
    for (i = 0; i < NS; i++) {
        
        // dummy elemental compositions
        dumC = randu();
        dumN = randu();
        dumP = randu();
        
        sumdum = dumC + dumN + dumP;
        
        // Assigning elemental compositions
        percentC[i] = dumC / sumdum;
        percentN[i] = dumN / sumdum;
        percentP[i] = dumP / sumdum;
        
        // calculate euclidian distance form the center of 1/3 1/3 1/3
        Cent[i] = EucDistFunction((1.0/3.0), (1.0/3.0), (1.0/3.0), percentC[i], percentN[i], percentP[i]);
        //printf("%f\t%f\t%f\t%f\t%f\n", Cent[i], (1.0/3.0), percentC[i], percentN[i], percentP[i]);
        // randomly assign primary producer of consumer
        //prim[i] = gsl_ran_binomial(r,PROPPRIMS,1);
        
        if ( prim[i] )
        numPrims++;
        else
        numCons++;
    }
}


void InitializationFunction(unsigned long int SeedSet)
{
    
    int i;
    
    // initialize random number generator
    gsl_rng_env_setup(); // set up the environment variables for the RNG
    T = gsl_rng_mt19937; // Mersenne Twister// for default you can say T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    gsl_rng_set(r, SeedSet);
    
    // Assigning Trophic Groups
    /*numPrimsForced = (int) (PROPPRIMS * (double) NS + 0.5);
     numHerbsForced = (int) (PROPHERBS * (double) (NS - numPrimsForced) + 0.5);
     numPredsForced = (int) NS - numPrimsForced - numHerbsForced;*/
    
    numPrimsForced = (int) (PROPPRIMS * (double) NS + 0.5);
    numHerbsForced = (int) (PROPHERBS * (double) NS + 0.5);
    numPredsForced = (int) (PROPPREDS * (double) NS + 0.5);
    
    //printf("%i,%i,%i\n", numPrimsForced, numHerbsForced, numPredsForced);
    //printf("%i,%i,%i\n", numPrimsForced, numHerbsForced+numPrimsForced, numPredsForced+numPrimsForced+numHerbsForced);
    
    for (i = 0; i < NS; i++) {
        if (i < numPrimsForced) {
            
            TrophicLevel[i] = 0;
            prim[i] = 1;
            pred[i] = 0;
        }
        else if (i < numHerbsForced + numPrimsForced) {
            
            TrophicLevel[i] = 1;
            prim[i] = 0;
            pred[i] = 0;
        }
        else {
            
            TrophicLevel[i] = 2;
            prim[i] = 0;
            pred[i] = 1;
        }
        
        //printf("%i %i\n", i, TrophicLevel[i]);
    }
    
    /*for (i = 0; i < NS; i++) {
     if (i <= numPrimsForced) {
     
     TrophicLevel = 1;
     }
     
     if (i > numPrimsForced && i <= numHerbsForced) {
     
     TrophicLevel = 2;
     }
     
     if (i > numHerbsForced) {
     
     TrophicLevel = 3;
     }
     
     printf("%i\n", TrophicLevel);
     }
     */
    
    for (i = 0; i < NS; i++) {
        
        
        if (prim[i] == 0) {
            // Initializing popuation sizes
            do {
                
                // population size from normal distribution
                PopSize[i] = MEANINITConsPOPSIZE + gsl_ran_gaussian (r, SDINITConsPOPSIZE);
                
                // if popsize is less than zero, choose again
            } while (PopSize[i] < 0);
        }
        
        if (prim[i] == 1) {
            // Initializing popuation sizes
            do {
                
                // population size from normal distribution
                PopSize[i] = MEANINITPrimPOPSIZE + gsl_ran_gaussian (r, SDINITPrimPOPSIZE);
                K[i] = Kcap + gsl_ran_gaussian (r, SDINITPrimPOPSIZE); // changeing (r,0) to (r,SDInitPrimPOPSIZE)
                
                // if popsize is less than zero, choose again
            } while (PopSize[i] < 0);
        }
    }
    
    // print to file every printTime time steps
    if ( t == 0 ) {
        
        // loop throught to print after process is over
        for (i = 0; i < NS; i++ ) {
        
            //printf("InitF %i %f %f %f %f %f %f \n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i]);
            
            fprintf(InitFile, "%i,%f,%f,%f,%f,%f,%f,%i,%i\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], prim[i],TrophicLevel[i]);
            
            
            if (prim[i] == 0) {
                
                fprintf(ConsInitFile, "%i,%f,%f,%f,%f,%f,%f,%i,%i\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], prim[i],TrophicLevel[i]);
            }
            if (prim[i] == 1) {
                
                fprintf(PrimInitFile, "%i,%f,%f,%f,%f,%f,%f,%i,%i\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], prim[i],TrophicLevel[i]);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

// Choosing three Prey items Fuunction

void chooseFunction(int nChoices) {
    int tempArray[nChoices];
    int focali, i, indexChoices[NS];
    int rechoose;
    int choiceArray[(NS-1)];
    int count, startChoice, nFoodOptions;
    
    TriangleFileInit = fopen("TriangleFileInit.csv","w");
    
    // using gsl
    // step 1 build the choice array
    
    for ( i = 0; i < NS; i++ )
    indexChoices[i] = i;
    
    for ( i = numPrimsForced; i < NS; i++ ) {
        if ( TrophicLevel[i] == 1 ) { // herbivores
            startChoice = 0;
            nFoodOptions = numPrimsForced;
        }
        else if ( TrophicLevel[i] == 2 ) { // preds
            startChoice = numPrimsForced;
            nFoodOptions = numHerbsForced;
        }
        else {
            fprintf(stdout, "\nerror in ChooseFunction(): undefined trophic level:\n");
            fprintf(stdout, "\tTrophicLevel[%i] = %i\n", i, TrophicLevel[i]);
            exit(-1);
        }
        
        gsl_ran_choose(r, &tempArray[0], nChoices, &indexChoices[startChoice], nFoodOptions, sizeof(int));
        
        aPrey[i] = tempArray[0];
        bPrey[i] = tempArray[1];
        cPrey[i] = tempArray[2];
        
        //printf("%i, %i, %i, %i, %i\n", i, TrophicLevel[i], aPrey[i], bPrey[i], cPrey[i]);
        
    }
    //exit(0);
    
    
    //    for (focali=0; focali < NS; focali++) {
    //
    //        count = 0;
    //        for ( i=0; i < NS; i++ ) {
    //            if ( focali != i ) {
    //                choiceArray[count] = i;
    //                count++;
    //            }
    //        }
    //        //check: unit testing building the choice
    //        for ( i=0; i < (NS-1); i++ ) {
    //            if ( focali == choiceArray[i] || choiceArray[0]==choiceArray[1] || choiceArray[0]==choiceArray[2] || choiceArray[1]==choiceArray[2] ) {
    //                printf("Bad algorithm!! Exiting with status -1\n\n");
    //                exit(-1);
    //            }
    //        }
    //
    //        // step 2: use gsl_ran_choose
    //        if (TrophicLevel[focali] == 2) {
    //
    //            gsl_ran_choose(r, &tempArray[0], nChoices, &choiceArray[0], (numPrimsForced-1), sizeof(int));
    //        }
    //
    //        if (TrophicLevel[focali] == 3) {
    //
    //            gsl_ran_choose(r, &tempArray[0], nChoices, &choiceArray[numHerbsForced], ((NS-1)-numPrimsForced), sizeof(int));
    //        }
    //
    //        /*
    //        else {
    //            gsl_ran_choose(r, &tempArray[0], nChoices, &choiceArray[0], (NS-1), sizeof(int));
    //        }*/
    //
    //        //gsl_ran_choose(r, &tempArray[0], nChoices, &choiceArray[0], (NS-1), sizeof(int));
    //
    //        // grab prey items from temp array
    //
    //        aPrey[focali] = tempArray[0];
    //        bPrey[focali] = tempArray[1];
    //        cPrey[focali] = tempArray[2];
    //    }
    for ( i=0; i<NS; i++ ) {
        
        float p0x[NS];
        float p1x[NS];
        float p2x[NS];
        
        float p0y[NS];
        float p1y[NS];
        float p2y[NS];
        int Triangle[NS];
        
        float Area[NS];
        float s[NS];
        float tt[NS];
        float p[NS];
        
        float px[NS];
        float py[NS];
        
        px[i] = percentC[i];
        py[i] = percentN[i];
        
        // percent carbon composition of prey
        percentC_aPrey[i] = percentC[aPrey[i]];
        percentC_bPrey[i] = percentC[bPrey[i]];
        percentC_cPrey[i] = percentC[cPrey[i]];
        // for traingel algorithm
        p0x[i] = percentC[aPrey[i]];
        p1x[i] = percentC[bPrey[i]];
        p2x[i] = percentC[cPrey[i]];
        
        // percent nitrogen composition of prey
        percentN_aPrey[i] = percentN[aPrey[i]];
        percentN_bPrey[i] = percentN[bPrey[i]];
        percentN_cPrey[i] = percentN[cPrey[i]];
        // for triangle algorithm
        p0y[i] = percentN[aPrey[i]];
        p1y[i] = percentN[bPrey[i]];
        p2y[i] = percentN[cPrey[i]];
        
        // percent phosphorus composition of prey
        percentP_aPrey[i] = percentP[aPrey[i]];
        percentP_bPrey[i] = percentP[bPrey[i]];
        percentP_cPrey[i] = percentP[cPrey[i]];
        
        // determines if the composition is in or outside the triangle
        Area[i] = 0.5 * ( (-p1y[i] * p2x[i]) + (p0y[i] * (-p1x[i] + p2x[i] )) + (p0x[i] * (p1y[i]  - p2y[i] )) + (p1x[i] * p2y[i]) );
        //printf("\narea = %f, %f\n", Area, (p0y * (-p1x  + p2x )));
        
        s[i]  = (1.0 / (2.0 * Area[i])) * (p0y[i] * p2x[i]  - p0x[i] * p2y[i] + (p2y[i]  - p0y[i] ) * px[i]  + (p0x[i]  - p2x[i] ) * py[i] );
        tt[i]  = (1.0 / (2.0 * Area[i])) * (p0x[i] * p1y[i]  - p0y[i] * p1x[i] + (p0y[i]  - p1y[i] ) * px[i]  + (p1x[i]  - p0x[i] ) * py[i] );
        
        if (s[i]  >= 0 && s[i]  <= 1 && tt[i] >= 0 && tt[i]  <= 1 && (s[i]  + tt[i] ) <= 1) {
            
            Triangle[i]  = 1;
        }
        
        else {
            
            Triangle[i]  = 0;
        }
        
        fprintf(TriangleFileInit,"%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i,%f,%f\n", i, Triangle[i],s[i],tt[i],p0x[i],p1x[i],p2x[i],p0y[i],p1y[i],p2y[i],aPrey[i],bPrey[i],cPrey[i],Cent[i],PopSize[i]);
    }
    // now you have 3 random choices stored in the tempArray
    /*for ( count = 0; count < 3; count++ )
     printf("%i %i \n", count, tempArray[count]);
     */
    
    /*for ( i = 0; i < NS; i++ )
     printf("ChooseF %i %i %i %i \n", i, aPrey[i], bPrey[i], cPrey[i]); */
    //printf("afrer print in choose \n");

    fclose(TriangleFileInit);

}

//////////////////////////////////////////////////////////////////////////////////////////

// funciton to calculate euclidian distance between two vectors
double EucDistFunction(double x1, double y1, double z1, double x2, double y2, double z2) {
    
    // returns a float when function is called
    return sqrt(((x1-x2) * (x1-x2)) + ((y1-y2) * (y1-y2)) + ((z1-z2) * (z1-z2)));
    
}

//////////////////////////////////////////////////////////////////////////////////////////


void OptForagingFunction(void){
    
    //printf("entered OF \n");
    int i;
    double effa;
    double effb;
    double effc;
    double effMin = 0.0;
    double effMax = 1.0;
    double effMaxa, effMaxb;
    //printf("all declared \n");
    
    for (i=0; i < NS; i++) {
        
        /*printf("%i %f %f %f %f %i \n", i, percentC_aPrey[i], percentN_aPrey[i], percentP_aPrey[i], (percentC_aPrey[i] + percentN_aPrey[i] + percentP_aPrey[i]),prim[i]);*/
        
        // if it is a primary producer
        if (prim[i] == 1) {
            
            // primary producers have no foraging effort
            OptEffa[i] = 0.0;
            OptEffb[i] = 0.0;
            OptEffc[i] = 0.0;
            // set precent element to zero
            percentC_Diet[i] =0.0;
            percentN_Diet[i] = 0.0;
            percentP_Diet[i] = 0.0;
            // set diet composition to zero
            percentC_Diet_Norm[i] = 0.0;
            percentN_Diet_Norm[i] = 0.0;
            percentP_Diet_Norm[i] = 0.0;
            EdistOpt[i] = 0.0;
            Edist[i] = 0.0;
            
        }
        // esle do this for consumers
        else if (PopSize[i] > 2.0) { // changed from >0 becasue dead ones Edist=1 // else if (PopSize[i] >= 0.0)
            
            EdistOpt[i] = 1000000000.0;  // changed from really high number
            //Optrawgains[i] = 100000000.0;
            
            OptEffa[i] = 0.0;
            OptEffb[i] = 0.0;
            OptEffc[i] = 0.0;
            
            // calculate optimal foraging efforts for all consumers
            if ( PopSize[aPrey[i]] <= 0.0 )
            effMaxa = 0.0;
            
            else
            effMaxa = effMax;
            effa = 0.0;
            effb = 0.0;
            effc = 0.0;
            for (effa=effMin; effa<=effMaxa; effa+=EFFORTINCR) {
                if ( PopSize[bPrey[i]] <= 0.0 )
                effMaxb = 0.0;
                else
                effMaxb = (1.0 - effa);
                
                for (effb=effMin; effb<=effMaxb; effb+=EFFORTINCR) {
                    if ( PopSize[cPrey[i]] <= 0.0 )
                    effc = 0.0;
                    else
                    effc = 1.0 - (effa + effb);
                    
                    //    if (effc < 0.00001) {
                    //    //normalize efforts if problem
                    //        effc = 0.0;
                    //        effa = effa/(effa + effb);
                    //        effb = effb/(effa + effb);
                    //
                    //                    }
                    // redo eff if above 1
                    if (effc > 1.0) effc = 1.0;
                    
                    // percent carbon in diet
                    percentC_Diet[i] = effa * PopSize[aPrey[i]] * percentC_aPrey[i] + effb * PopSize[bPrey[i]] * percentC_bPrey[i] + effc * PopSize[cPrey[i]] * percentC_cPrey[i];
                    // percent nitrogen in diet
                    percentN_Diet[i] = effa * PopSize[aPrey[i]] * percentN_aPrey[i] + effb * PopSize[bPrey[i]] * percentN_bPrey[i] + effc * PopSize[cPrey[i]] * percentN_cPrey[i];
                    // percent phosphorus in diet
                    percentP_Diet[i] = effa * PopSize[aPrey[i]] * percentP_aPrey[i] + effb * PopSize[bPrey[i]] * percentP_bPrey[i] + effc * PopSize[cPrey[i]] * percentP_cPrey[i];
                    
                    // if no foraging happens, leave populqtion size the same
                    if (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i] == 0) {
                        
                        PopSize[i] = PopSize[i];
                        //EdistOpt[i] = NAN;
                        
                    }
                    else {
                        // normalize diets so they sum to 1
                        percentC_Diet_Norm[i] = percentC_Diet[i] / (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]);
                        
                        // normalize diets so they sum to 1
                        percentN_Diet_Norm[i] = percentN_Diet[i] / (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]);
                        
                        // normalize diets so they sum to 1
                        percentP_Diet_Norm[i] = percentP_Diet[i] / (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]);
                        
                        // Calculate Euclidian Distance between diet and composition
                        Edist[i] = EucDistFunction(percentC_Diet_Norm[i], percentN_Diet_Norm[i], percentP_Diet_Norm[i], percentC[i], percentN[i], percentP[i]);
                        
                        
                        //rawgains[i] = ((Gains[i] * exp (-THETA * (Edist[i]*Edist[i]))) - TotalKilled[i] - BASALMR * PopSize[i]);
                        
                        // update optimal edist if current edist is better
                        //if (rawgains[i] < Optrawgains[i]) {
                        if (Edist[i] < EdistOpt[i]) {
                            //Optrawgains[i] = rawgains[i];
                            EdistOpt[i] = Edist[i];
                            OptEffa[i] = effa;
                            OptEffb[i] = effb;
                            OptEffc[i] = effc;
                            
                        }
                    } // end of else statement
                }
            }
            if ( (OptEffa[i] > 0.0 && PopSize[aPrey[i]] <= 0.0) || (OptEffb[i] > 0.0 && PopSize[bPrey[i]] <= 0.0) || (OptEffc[i] > 0.0 && PopSize[cPrey[i]] <= 0.0) ) {
                //printf("\nError in optimal foraging: species %i, \naPrey = %i, bPrey = %i, cPrey = %i\nPopSize[aPrey[i]] = %f, PopSize[bPrey[i]] = %f, PopSize[cPrey[i]] = %f\nPopSize[i] = %f\nOptEffa[i] = %f, OptEffb[i] = %f, OptEffc[i] = %f\n", i, aPrey[i], bPrey[i], cPrey[i], PopSize[aPrey[i]], PopSize[bPrey[i]], PopSize[cPrey[i]], PopSize[i], OptEffa[i], OptEffb[i], OptEffc[i]);
                exit(-1);
            }
        } // end of else statement
        
        //printf("OptForgF %i %f %f %f %f %i %f %f \n", i, Edist[i], (percentC_Diet_Norm[i] + percentN_Diet_Norm[i] + percentP_Diet_Norm[i]), (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]), (OptEffa[i] + OptEffb[i] + OptEffc[i]), prim[i],Edist[i],EdistOpt[i]);
        // print to file every printTime time steps
        if ( t == (T_TOTAL - 1) ) {
            
            fprintf(OptForgFile,"%i,%f,%f,%f,%f,%i,%f,%f,%f,%f\n", i, Edist[i], (percentC_Diet_Norm[i] + percentN_Diet_Norm[i] + percentP_Diet_Norm[i]), (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]), (OptEffa[i] + OptEffb[i] + OptEffc[i]), prim[i],EdistOpt[i],OptEffa[i],OptEffb[i],OptEffc[i]);
            
            if (prim[i] == 0) {
                
                fprintf(ConsOptForgFile,"%i,%f,%f,%f,%f,%i,%f,%f,%f,%f\n", i, Edist[i], (percentC_Diet_Norm[i] + percentN_Diet_Norm[i] + percentP_Diet_Norm[i]), (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]), (OptEffa[i] + OptEffb[i] + OptEffc[i]), prim[i],EdistOpt[i],OptEffa[i],OptEffb[i],OptEffc[i]);
            }
            
            if (prim[i] == 1) {
                
                fprintf(PrimOptForgFile,"%i,%f,%f,%f,%f,%i,%f\n", i, Edist[i], (percentC_Diet_Norm[i] + percentN_Diet_Norm[i] + percentP_Diet_Norm[i]), (percentC_Diet[i] + percentN_Diet[i] + percentP_Diet[i]), (OptEffa[i] + OptEffb[i] + OptEffc[i]), prim[i],EdistOpt[i]);
            }
            
        }
        
    } // end in i loop
    
} // end of function

//////////////////////////////////////////////////////////////////////////////////////////




//GainLossFunction calculating raw gains and losses

void GainLossFunction(void) {
    
    int i, j, l, m;
    
    if (t < (T_TOTAL + 1) /** StartCollectTime*/) {
        
        for (i=0; i<NS; i++) {
            for (j=0; j<NS; j++) {
                
                GainLossAll[i][j] = 0.0;
            }
        }
    }
    
    // for all rows i, do each column j
    for (i=0; i<NS; i++) {
        
        Gains[i] = 0.0;
        TotalKilled[i] = 0.0;
        
        for(j=0; j<NS; j++) {
            
            // Reset GainLoss to zero everytime
            GainLoss[i][j] = 0.0;
            EffortMatrix[i][j] = 0.0;
            
            if ( !prim[j] ) {
                
                // Change GainLoss for prey items a, b, and c
                if (aPrey[j] == i) {
                    
                    GainLoss[i][j] = deltaT * AttackRate * OptEffa[j] * PopSize[j] * PopSize[i];
                    EffortMatrix[i][j] = OptEffa[j];
                    
                    GainLossAll[i][j] = GainLossAll[i][j] + GainLoss[i][j];
                }
                // Change GainLoss for prey items a, b, and c
                else if (bPrey[j] == i) {
                    
                    GainLoss[i][j] = deltaT * AttackRate * OptEffb[j] * PopSize[j] * PopSize[i];
                    EffortMatrix[i][j] = OptEffb[j];
                    // Sum of interatcions through time
                    GainLossAll[i][j] = GainLossAll[i][j] + GainLoss[i][j];
                }
                // Change GainLoss for prey items a, b, and c
                else if (cPrey[j] == i) {
                    
                    GainLoss[i][j] = deltaT * AttackRate * OptEffc[j] * PopSize[j] * PopSize[i];
                    EffortMatrix[i][j] = OptEffc[j];
                    // Sum of interatcions through time
                    GainLossAll[i][j] = GainLossAll[i][j] + GainLoss[i][j];
                }
                // Stopping species from eating nonexistnant prey itmes
                if (GainLoss[i][j] > PopSize[i]) {
                    // Sum of interatcions through time
                    GainLoss[i][j] = PopSize[i];  // added in the scaler here???
                }
            }
            /*if ( check && GainLoss[i][j] < 0.0001 && !prim[j] ) {
             printf("%f %f %f %f %f %f %i %i %i %i %i %i %i", GainLoss[i][j], OptEffa[j], OptEffb[j], OptEffc[j], PopSize[j], PopSize[i], i, j, aPrey[j], bPrey[j], cPrey[j],prim[j],prim[i]); // printing GainLoss all rows
             exit(0);
             }*/
            // print to file every printTime time steps
            if ( t == (T_TOTAL - 1) ) {
                
                fprintf(GainLossFile,"%f,", GainLoss[i][j]); // printing GainLoss all rows
                fprintf(GainLossAllFile,"%f,", GainLossAll[i][j]); // printing GainLoss all rows
                fprintf(EffortMatrixFile,"%f,", EffortMatrix[i][j]);
                
            }
        }
        //printf("\n"); // then carriage return
        // print to file every printTime time steps
        if ( t == (T_TOTAL - 1) ) {
            
            fprintf(GainLossFile,"\n"); // then carriage return
            fprintf(GainLossAllFile,"\n"); // then carriage return
            fprintf(EffortMatrixFile,"\n");
            
        }
        
        for (j=0; j < NS; j++) {
            
            // summing all predation on each species
            TotalKilled[i] = TotalKilled[i] + GainLoss[i][j];
            
        }
        
        if (TotalKilled[i] > PopSize[i]) {
            
            for (j=0; j<NS; j++) {
                
                // correcting all the predation that exceeds popsizes
                GainLoss[i][j] =  PopSize[i] * (GainLoss[i][j] / TotalKilled[i]);
                
            }
            
            TotalKilled[i] = PopSize[i];
            
        } // calculating foraging gains for each species
        
        
        for (m=0; m<NS; m++) {
            
            Gains[i] = Gains[i] + GainLoss[m][i];
            
        } // end of m loop
        
    } // end of i loop
    
    // calc the sum of gains for consumers & if j
    // print out the predation and foraging gains
    /*for (i=0; i<NS; i++) {
     
     printf("TotalKilled[i] %i %f %f \n", i, TotalKilled[i], Gains[i]);
     }*/
}

//////////////////////////////////////////////////////////////////////////////////////////


// Population dynamics function
void PopulationUpdateFunction(void) {
    
    
    double PopSizeTemp;
    int i;
    
    Stable_Trigger = 0;
    
    // for each species i
    for (i=0; i<NS; i++) {
        
        // set up temp population size
        PopSizeTemp = PopSize[i];
        
        if (prim[i] == 0) {
            
            // population dynamics for consumers
            PopSizePrevious[i] = PopSize[i];
            
            //Original matcing function
            PopSize[i] = PopSize[i] + (( GAINSSCALER * Gains[i] * exp (-THETA * (Edist[i]*Edist[i]))) - TotalKilled[i] - deltaT * BASALMR * PopSize[i]); // multiplied BASALMR by deltaT
            
            //Modified version of matching function
            //PopSize[i] = PopSize[i] + (( GAINSSCALER * Gains[i] * exp (-THETA * (Edist[i]))) - TotalKilled[i] - deltaT * BASALMR * PopSize[i]); // multiplied BASALMR by deltaT
            
        }
        else {
            
            // population dynamics for primary producers
            PopSizePrevious[i] = PopSize[i];
            
            PopSize[i] = PopSize[i] + deltaT * ((rPP * PopSize[i] * ((K[i] - PopSize[i])/K[i])) - TotalKilled[i]);
        }
        
        PopSizeChange[i] = fabs((PopSizePrevious[i] - PopSize[i]) / PopSizePrevious[i]);
        
        TotalPopSizeChange = TotalPopSizeChange + PopSizeChange[i];
        
        if (t < MinTimeStable || (PopSize[i] > (MeanPopSize/10) && PopSizeChange[i] > ChangeCritical) ) {
            
            Stable_Trigger = 1;
        }
        
        else {
            
            for (i=0; i < NS; i++) {
                
                fprintf(LastTimeFile, "%i,%i,%f,%f\n", t,i,PopSizeChange[i],(PopSizePrevious[i] - PopSize[i]));
                
            }
        }
        
    }
    
    // run through population sizes and stop any negaitives
    for (i=0; i<NS; i++) {
        
        
        // if really low population size, let them go extinct
        if (PopSize[i] < 2.0) {
            
            PopSize[i] = 0.0;
            
            
            //fprintf(ExtinctionLogFile, "%i,%i,%f,%i,%f,%i,%f,%f,%f,%i,%f\n", i, aPrey[i],PopSize[aPrey[i]],bPrey[i],PopSize[bPrey[i]],cPrey[i],PopSize[cPrey[i]], Cent[i], PopSize[i], t, Edist[i]);
            
            
            if (PopSizePrevious[i] > 0.0 && PopSize[i] == 0.0) {
                
                fprintf(ExtinctionLogFile, "%i,%i,%f,%f,%f,%i,%f,%f,%i,%f,%f,%i,%f,%f,%f,%f,%i,%f\n", TrophicLevel[i],i,Gains[i],TotalKilled[i],((Gains[i]*GAINSSCALER)-TotalKilled[i]),aPrey[i],OptEffa[i],PopSize[aPrey[i]],bPrey[i],OptEffb[i],PopSize[bPrey[i]],cPrey[i],OptEffc[i],PopSize[cPrey[i]], Cent[i], PopSize[i], t, Edist[i]);
                
            }
        }
        // print to file every printTime time steps
        if ( t > (MinTimeStable -1) ) {
            
            fprintf(PopUpdateFile, "%i,%f,%f,%f,%f,%f,%f,%i,%f\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], t, Edist[i]);
            
            if (prim[i] == 0) {
                
                fprintf(ConsPopUpdateFile, "%i,%f,%f,%f,%f,%f,%f,%i,%f\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], t, Edist[i]);
            }
            if (prim[i] == 1) {
                
                fprintf(PrimPopUpdateFile, "%i,%f,%f,%f,%f,%f,%f,%i,%f\n", i, percentC[i], percentN[i], percentP[i], (percentC[i] + percentN[i] + percentP[i]), Cent[i], PopSize[i], t, Edist[i]);
            }
        }
    }
    
    /*if (Stable_Trigger == 0){
     
     t=T_TOTAL;
     } */
    
    
} // end of function

//////////////////////////////////////////////////////////////////////////////////////////


// new function usage() added by SMF 10.14.14
void usage(char *s)
{
    fprintf(stderr,  "\nUsage: %s [options]\n\n", s);
    
    fprintf(stderr,  "\tNOTE: There is no automated error checking on user inputs.\n");
    fprintf(stderr,  "\tProgram behavior is unpredictable (not defined) for\n");
    fprintf(stderr,  "\tinputs that do not conform to guidelines below.\n\n");
    
    fprintf(stderr,  "\tNOTE: There is also no automated error checking on this\n");
    fprintf(stderr,  "\toutput.  This usage() function was last updated 10.14.14.\n");
    fprintf(stderr,  "\tHence, if the code was subsequently updated, the descriptions\n\tbelow may be out of date.\n\n");
    
    
    fprintf(stderr,  "\t-d <step size>\n\t\tSize of discrete time step, deltaT\n");
    fprintf(stderr,  "\t\tDefault is %f\n\n", deltaT_DEFAULT);
    
    fprintf(stderr,  "\t-K <value>\n\t\tCarrying capacity for primary producers, K\n");
    fprintf(stderr,  "\t\tDefault is %f\n\n", K_DEFAULT);
    
    fprintf(stderr,  "\t-T <integer>\n\t\tTotal number of discrete time steps, T_TOTAL, K\n");
    //fprintf(stderr,  "\t\tDefault is %i\n\n", T_TOTAL_DEFAULT);
}


///////////////////////////////////////////////////////////////////////////////////////


// Calculate trophic level using shortest chain method
void calcChainLength(void) {
    
    int preyChainLength[NS];
    int count;
    int i, j;
    int looking = 1, focalLevel;
    int deadOnes;
    int totNum;
    int numAtLevel;
    int PredNoPrey = 0;
    
    for ( i = 0; i < NS; i++ ) {
        if ( prim[i] && PopSize[i] > 0.0)
        chainLength[i] = 0;
        else
        chainLength[i] = -1;
    }
    
    focalLevel = 1;
    count = 1;
    
    while ( looking ) {
        for ( i = 0; i < NS; i++ ) {
            
            if ( chainLength[i] == -1 && PopSize[i] > 0.0 ) {
                // figure out WHO i eats ONLY if you have not yet found i's chain length
                // only count prey items with pop abundance > 0.0
                
                if (PopSize[aPrey[i]] * OptEffa[i] + PopSize[bPrey[i]] * OptEffb[i] + PopSize[cPrey[i]] * OptEffc[i] > 0.0) {
                    preyChainLength[0] = chainLength[aPrey[i]];
                    preyChainLength[1] = chainLength[bPrey[i]];
                    preyChainLength[2] = chainLength[cPrey[i]];
                    
                    // check its food items to see if chainLength[] == (focalLevel - 1)
                    for (j=0; j<3; j++) {
                        if ( preyChainLength[j] == (focalLevel-1) ) {
                            chainLength[i] = focalLevel;
                        }
                    }
                }
                else {
                    (chainLength[i] = -2);
                    PredNoPrey++;
                }
            }
        } // you have checked all for focal levels <= focalLevel
        // NOW see if any remain to check;
        looking = 0;
        for ( i = 0; i < NS; i++ ) {
            if ( chainLength[i] == -1 && PopSize[i] > 0.0 ) {
                looking = 1;
            }
        }
        focalLevel++;
        count++;
        
        if ( count > (NS + 1) ) {
            printf("\nError in calcChainLength(): infinite loop!\n");
            printf("%i   %i    %i %f %i %i %i\n", i, chainLength[i], prim[i], PopSize[i],aPrey[i],bPrey[i],cPrey[i]);
            exit(-1);
        }
    }
    
    for ( i = 0; i < NS; i++ ) {
        printf("%i %i %i %f %i %i %i\n", i, chainLength[i], prim[i], PopSize[i],aPrey[i],bPrey[i],cPrey[i]);
        fprintf(TrophicFile,"%i,%i,%i,%f,%i,%i,%i\n", i,chainLength[i],prim[i],PopSize[i],aPrey[i],bPrey[i],cPrey[i]);
        if ( PopSize[i] > 0.0 && chainLength[i] < 0 && chainLength[i] > -2 ) {
            printf("PROBLEM!!! PopSize > 0 and chainLength < 0 \n");
        }
    }
    
    totNum = 0;
    for ( i = 0; i < focalLevel; i++ ) {
        numAtLevel = 0;
        for ( j = 0; j < NS; j++ ) {
            if ( chainLength[j] == i )
            numAtLevel++;
        }
        printf("%i species at level %i\n", numAtLevel, i);
        totNum = totNum + numAtLevel;
    }
    
    deadOnes = 0;
    for ( i = 0; i < NS; i++ ) {
        if ( !(PopSize[i] > 0.0) )
        deadOnes++;
    }
    
    if ( (totNum + deadOnes + PredNoPrey) != NS ) {
        printf("ERROR!!! Does not sum to NS, %i \n", totNum+deadOnes+PredNoPrey);
    }
}

// Print parameters function
void PrintParams(void) {
    
    fprintf(ParamFile,"%f,%i,%i,%f,%f,%i,%i,%i,%i,%f,%f,%f,%i,%i,%f,%f,%f,%i,%i,%f,%f,%f\n", deltaT,T_TOTAL,MEANINITPrimPOPSIZE,rPP,PROPPRIMS,MEANINITPrimPOPSIZE,SDINITPrimPOPSIZE,MEANINITConsPOPSIZE,SDINITConsPOPSIZE,BASALMR,EFFORTINCR,THETA,ADDEDSPECIES,REMOVEDSPECIES,GAINSSCALER,LOWSPACE,HIGHSPACE,PARAMSTEPS,PRINTTIME,MeanPopSize,MeanPrimPopSize,MeanConsPopSize);
}

// Function to drop the most elementally extreme species
void DropExtremeSpecies(void) {
    
    int ppCurrentExtremeSpecies;
    int ppExtremeSpecies;
    double ppCentMax = 0;
    double ppCentMin = 100;
    int i;
    
    for (i = 0; i < NS; i++)
    
    if (prim[i] == 1 && PopSize[i] > 10 && TotalKilled[i] > 0 && Cent[i] > ppCentMax) {
        
        ppCentMax = Cent[i];
        ppCurrentExtremeSpecies = i;
    }
    
    ppExtremeSpecies = ppCurrentExtremeSpecies;
    
    PopSize[ppExtremeSpecies] = 0;
    
}
