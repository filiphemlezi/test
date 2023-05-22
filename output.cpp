#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
    #include <stdlib.h>
        #include <thread>
            #include <random>
                #include <mutex>
                    #include <cstring>

                        #include "output.hpp"

const float V = ( 2.00537614e-6 ) ;
const float dt = ( 50.0 ) ;
const float K0 = ( 0.0000444 ) ;
const float m0 = ( 1.67261e-27 ) ;
const float q = ( 1.60219e-19 ) ;
const float c = ( 2.99793e8 ) ;
const float T0 = ( ( 1.67261e-27 ) * ( ( 2.99793e8 ) * ( ( 2.99793e8 ) / ( ( 1.60219e-19 ) * ( 1e9 ) ) ) ) ) ;

                        float getTkinInjection(unsigned long long id, float from, float to, int numberOfBins)
                        {
                        float step = (to - from) / numberOfBins;
                        return from + ((id % numberOfBins) * step);
                        }

                        float getSolarPropInjection(unsigned long long state, int numberOfBins)
                        {
                        int modulo = state % numberOfBins;
                        return 0.01f * powf((1.0f + 0.5f), modulo);
                        }


                        void preSimulation(preSimulationStructure *output, int padding)
                        {
                        int id;
                        int m, mm;
                        for (m = 0; m < 101; m++)
                        {
                        for (mm = 0; mm < 250; mm++)
                        {
                        id = padding+mm;
            float Tkin = ( getTkinInjection(id,0.0001,101,10000) ) ;
            float rigSqrt = ( ( Tkin ) * ( ( Tkin ) + ( ( 2.0 ) * ( T0 ) ) ) ) ;
            float Rig = ( sqrtf(rigSqrt) ) ;
            float p = ( ( Rig ) * ( ( 1e9 ) * ( ( q ) / ( c ) ) ) ) ;
            output[id].Tkininj = ( Tkin ) ;
            output[id].pinj = ( p ) ;
                        }
                        }
                        }

                        void simulation(preSimulationStructure *input, simulationStructure *output, int padding)
                        {
                        int id;
                        thread_local std::random_device rd{};
                        thread_local std::mt19937 generator(rd());
                        thread_local std::normal_distribution<float> distribution(0.0f, 1.0f);
                            int m, mm;
                            for (m = 0; m < 101; m++)
                            {
                            printf("%d\n",m);
                            for (mm = 0; mm < 250; mm++)
                            {
                            id = padding+mm;
            float r = ( 1.0 ) ;
            float p = ( input[id].pinj ) ;
            float beta ;
            float Rig ;
            float dr ;
            float pp ;
            float Tkin = ( input[id].Tkininj ) ;
            float TkinSqrt ;
            float betaSqrt ;
            float drSqrt ;
            while ( ( r )             < ( 100.0002 )) {
            betaSqrt = ( ( Tkin ) * ( ( Tkin ) + ( ( T0 ) + ( T0 ) ) ) ) ;
            beta = ( ( sqrtf(betaSqrt) ) / ( ( Tkin ) + ( T0 ) ) ) ;
            Rig = ( ( ( p ) * ( ( c ) / ( q ) ) ) / ( 1e9 ) ) ;
            pp = ( p ) ;
            p = ( ( p ) - ( ( 2.0 ) * ( ( V ) * ( ( p ) * ( ( dt ) / ( ( 3.0 ) * ( r ) ) ) ) ) ) ) ;
            drSqrt = ( ( 2.0 ) * ( ( K0 ) * ( ( beta ) * ( ( Rig ) * ( dt ) ) ) ) ) ;
            dr = ( ( ( ( V ) + ( ( 2.0 ) * ( ( K0 ) * ( ( beta ) * ( ( Rig ) / ( dt ) ) ) ) ) ) * ( dt ) ) + ( ( distribution(generator) ) * ( sqrtf(drSqrt) ) ) ) ;
            r = ( ( r ) + ( dr ) ) ;
            Rig = ( ( p ) * ( ( c ) / ( q ) ) ) ;
            TkinSqrt = ( ( ( T0 ) * ( ( T0 ) * ( ( q ) * ( ( q ) * ( ( 1e9 ) * ( 1e9 ) ) ) ) ) ) + ( ( q ) * ( ( q ) * ( ( Rig ) * ( Rig ) ) ) ) ) ;
            Tkin = ( ( ( sqrtf(TkinSqrt) ) - ( ( T0 ) * ( ( q ) * ( 1e9 ) ) ) ) / ( ( q ) * ( 1e9 ) ) ) ;
            Rig = ( ( Rig ) / ( 1e9 ) ) ;
            betaSqrt = ( ( Tkin ) * ( ( Tkin ) + ( ( T0 ) + ( T0 ) ) ) ) ;
            beta = ( ( sqrtf(betaSqrt) ) / ( ( Tkin ) + ( T0 ) ) ) ;
            if ( ( beta )             > ( 0.01 )) {
            if ( ( Tkin )             < ( 200.0 )) {
            if ( ( r )             > ( 100.0 )) {
            if ( ( ( r ) - ( dr ) )             < ( 100.0 )) {
            double w = ( ( ( m0 ) * ( ( m0 ) * ( ( c ) * ( ( c ) * ( ( c ) * ( c ) ) ) ) ) ) + ( ( p ) * ( ( p ) * ( ( c ) * ( c ) ) ) ) ) ;
            w = ( ( ( pow(w,-1.85) ) / ( p ) ) / ( 1e45 ) ) ;
            output[id].Tkininj = ( input[id].Tkininj ) ;
            output[id].Tkin = ( Tkin ) ;
            output[id].r = ( r ) ;
            output[id].w = ( w ) ;
            break;
            }
            }
            }
            }
            if ( ( beta )             < ( 0.01 )) {
            break;
            }
            if ( ( r )             < ( 0.3 )) {
            r = ( ( r ) - ( dr ) ) ;
            p = ( pp ) ;
            }
            }
                            }
                            }
                            printf("ID: %d done\n",(id/12)/250);
                            }

                            void postSimulation(simulationStructure *input, postSimulationStructure *output, int padding)
                            {
                            int id;
                            int m, mm;
                            for (m = 0; m < 101; m++)
                            {
                            for (mm = 0; mm < 250; mm++)
                            {
                            id = padding+mm;
            output[id].Tkininj = ( input[id].Tkininj ) ;
            output[id].Tkin = ( input[id].Tkin ) ;
            output[id].r = ( input[id].r ) ;
            output[id].w = ( input[id].w ) ;
                            }
                            }
                            }

                            void runModel()
                            {
                            unsigned int nthreads = std::thread::hardware_concurrency();
                            int N = 100000000;
                            int iterations = ceil((double)N/((double)nthreads*101*250));
                            printf("Spolu iteracii: %d\n",iterations);
                            preSimulationStructure* preSimulationStructureData, * preSimulationStructureLocal;
                            simulationStructure* simulationStructureData, * simulationStructureLocal;
                            postSimulationStructure* postSimulationStructureData, * postSimulationStructureLocal;
                            preSimulationStructureData = (preSimulationStructure*)malloc( sizeof(preSimulationStructure)*nthreads * iterations * 250);
                            preSimulationStructureLocal = (preSimulationStructure*)malloc(  sizeof(preSimulationStructure) *nthreads * iterations * 250);
                            simulationStructureData = (simulationStructure*)malloc( sizeof(simulationStructure) *nthreads * iterations * 250);
                            simulationStructureLocal = (simulationStructure*)malloc( sizeof(simulationStructure) *nthreads * iterations * 250);
                            postSimulationStructureData = (postSimulationStructure*)malloc( sizeof(postSimulationStructure)  *nthreads * iterations * 250);
                            postSimulationStructureLocal = (postSimulationStructure*)malloc( sizeof(postSimulationStructure) *nthreads * iterations * 250);
                            FILE* file = fopen("log.dat", "w");
                            std::vector<std::thread> threadsPreSimulation;
                                for (int k = 0; k < iterations; k++)
                                {
                                printf("Predsimulacna Iteracia: %d\n",k);
                                for (int i = 0; i < nthreads; i++)
                                {
                                threadsPreSimulation.emplace_back(std::thread(preSimulation, preSimulationStructureData,(nthreads*k+i)*250));
                                }
                                }
                                for (auto& th : threadsPreSimulation)
                                {
                                th.join();
                                }
                                std::vector<std::thread> threadsSimulation;
                                    for (int k = 0; k < iterations; k++)
                                    {
                                    printf("Simulacna Iteracia: %d\n",k);
                                    for (int i = 0; i < nthreads; i++)
                                    {
                                    threadsSimulation.emplace_back(std::thread(simulation, preSimulationStructureData, simulationStructureData, (nthreads*k+i)*250));
                                    }
                                    }
                                    for (auto& th : threadsSimulation) {
                                    th.join();
                                    }
                                    for (int k = 0; k < iterations; k++)
                                    {
                                    printf("Postsimulacna Iteracia: %d\n",k);
                                    std::vector<std::thread> threadsPostSimulation;
                                        for (int i = 0; i < nthreads; i++)
                                        {
                                        threadsPostSimulation.emplace_back(std::thread(postSimulation, simulationStructureData, postSimulationStructureData, (nthreads*k+i)*250));
                                        }
                                        for (auto& th : threadsPostSimulation) {
                                        th.join();
                                        }
                                        memcpy(postSimulationStructureLocal, postSimulationStructureData, sizeof(postSimulationStructure) * nthreads * iterations * 250);
                                        for (int i = (int)nthreads * k *250; i < (int)nthreads * (k+1) * 250; i++) {

                                        if (postSimulationStructureLocal[i].Tkininj != -1.0f) {
                                        fprintf(file, "%g,%g,%g,%g\n ", postSimulationStructureLocal[i].Tkininj, postSimulationStructureLocal[i].Tkin, postSimulationStructureLocal[i].r, postSimulationStructureLocal[i].w);
                                        }


                                        }
                                        }
                                        fclose(file);
                                        }


                                        int main() {
                                        runModel();
                                        return 0;
                                        }