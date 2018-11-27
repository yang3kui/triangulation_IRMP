#include "simulation.h"
#include <iostream>
#include <thread>
using namespace tri;
int main()
{
    Simulation* sim = nullptr;

    for (char i_sim = 'A'; i_sim <= 'D'; i_sim ++)
    {
        std::cout << "Configure " << i_sim << std::endl;
        switch(i_sim)
        {
        case 'A':
            sim = new SimulationConfig1();
            break;
        case 'B':
            sim = new SimulationConfig2();
            break;
        case 'C':
            sim = new SimulationConfig3();
            break;
        case 'D':
            sim = new SimulationConfig4();
            break;

        }

        sim->genSynData();
        sim->projection();
        sim->corrupt();
        sim->genBearingVec();
        for(int i = 0; i <= 2; i++)
        {
            // std::cout << "test: " << i << " begin!" << std::endl;
            sim->reconstruct(i);
            sim->computeErr();

            sim->print_result();
        }
    }
    std::cout << "done!" << std::endl;


    std::chrono::milliseconds time_span(1000);
    std::this_thread::sleep_for(time_span);
    delete sim;
    std::cin.get();
}
