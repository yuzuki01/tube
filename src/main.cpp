#include "tube.h"


int main(int argc, char **argv) {
    if (argc < 2) return 0;
    Tube tube(argv[1]);
    tube.initial();
    tube.output();

    while (true) {
        tube.do_step();
        if (tube.step % 10 == 0) {
            std::cout << "step: " << tube.step << "  solution time: " << tube.solution_time << std::endl;
            tube.output();
        }
        if (not tube.run_status()) {
            std::cout << "achieved stop time: " << tube.stop_time
                      << " step: " << tube.step << "  solution time: " << tube.solution_time << std::endl;
            break;
        }
    }
    tube.output();

    return 0;
}
