#include "tube.h"


int main(int argc, char **argv) {
    if (argc < 2) return 0;
    if (argc > 2) {
        std::vector<double> dp_Kn(argc - 2);
        for (int i = 2; i < argc; ++i) {
            dp_Kn[i - 2] = std::stod(argv[i]);
        }
        const int dp_num = int(dp_Kn.size());
        omp_set_num_threads(dp_num);
#pragma omp parallel for shared(dp_num, argv, std::cout, dp_Kn) default(none)
        for (int i = 0; i < dp_num; ++i) {
            Tube tube(argv[1], i + 1, dp_Kn[i]);
            tube.initial();
            tube.output();

            while (true) {
                tube.do_step();
                if (tube.step % 10 == 0) {
                    // std::cout << "step: " << tube.step << "  solution time: " << tube.solution_time << std::endl;
                    tube.output();
                }
                if (not tube.run_status()) {
#pragma omp critical
                    std::cout << "design point: " << i + 1 << " achieved stop time: " << tube.stop_time
                              << " step: " << tube.step << " solution time: " << tube.solution_time << std::endl;
                    break;
                }
            }
            tube.output();
        }
    } else {
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
    }

    return 0;
}
