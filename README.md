[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Penalty decomposition methods for second-best congestion pricing problems on large-scale networks

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Penalty decomposition methods for second-best congestion pricing problems on large-scale networks](https://doi.org/10.1287/ijoc.2023.0144) by Lei Guo, Wenxin Zhou, Xiaolei Wang, Hai Yang and Tijun Fan.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0144

https://doi.org/10.1287/ijoc.2023.0144.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Guo2024,
  author =        {Lei Guo and Wenxin Zhou and Xiaolei Wang and Hai Yang and Tijun Fan},
  publisher =     {INFORMS Journal on Computing},
  title =         {Penalty decomposition methods for second-best congestion pricing problems on large-scale networks},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0144.cd},
  url =           {https://github.com/INFORMSJoC/2023.0144},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0144},
}  
```

## Description

The second-best congestion pricing (SBCP) problem is one of the most challenging problems in transportation due to its  two-level hierarchical structure. In spite of various intriguing attempts for solving SBCP, existing solution methods are either heuristic without convergence guarantee or suitable for solving SBCP on small networks only. In this paper, we first reveal some convexity-based structural properties of the marginal value function reformation of SBCP and then, by effectively exploiting these structural properties, we propose two dedicated decomposition methods for solving SBCP on large-scale networks which are different from existing methods in that they avoid linearizing nonconvex functions. We establish the convergence of the two decomposition methods under commonly used conditions and provide the maximum number of iterations for deriving an approximate stationary solution. The computational experiments based on a collection of real road networks show that in comparison with three existing popular methods, the two proposed methods are capable of solving SBCP on larger-scale networks; and for instances that can be solved by existing methods, the two proposed methods are substantially faster.

This project contains two folders: `PD`, `RPD`.

- `PD`: This folder contains the data, source codes, makefile and results of the penalty decomposition methods.
- `RPD`: This folder contains the data, source codes, makefile and results of the relaxed penalty decomposition methods.

## Prerequisites

The codes are implemented under Ubuntu 18.04. Boost C++ 1.71.0 is also required for running the codes that is used in computing the shortest path tree.

## Building

In Linux/Ubuntu operating system, to build either `PD` or `RPD`, execute the following command in each folder.

```
make
```

To run the executable file, call the following command.

```
./main
```

Be sure to clean all the dependencies and executable files before building a different version of the code.

## Source File Description

- `PD\src\path_equilibration.h`: the header file for `PD\src\path_equilibration.cpp`.
- `PD\src\path_equilibration_rho.h`: the header file for `PD\src\path_equilibration_rho.cpp`.
- `PD\src\bbstep_PGM.h`: the header file for `PD\src\bbstep_PGM.cpp`.
- `PD\src\cost.hpp`: the header file containing the link performance function that used to solve either lower-level traffic assignment problem or Problem (11).
- `PD\src\path.hpp`: the header file containing the structure of the path.
- `PD\src\graph.hpp`: the header file containing the structure of the graph.
- `PD\src\io.hpp`: the header file to handle input/output streams.
- `PD\src\utils.hpp`: the header file containing auxiliary functions.
- `PD\src\path_equilibration.cpp`: the source codes for solving the lower-level traffic assignment problem with given $u$.
- `PD\src\path_equilibration_rho.cpp`: the source codes for solving Problem (11) (i.e., step 1 of Algorithm 1) to update $v$.
- `PD\src\bbstep_PGM.cpp`: the source codes for solving Problem (12) (i.e., step 2 of Algorithm 1) to update $u$.
- `PD\main.cpp`: the source codes of Algorithm PD.
- `RPD\src\path_equilibration.h`: the header file for `RPD\src\path_equilibration.cpp`.
- `RPD\src\path_equilibration_rho.h`: the header file for `RPD\src\path_equilibration_rho.cpp`.
- `RPD\src\bbstep_PGM.h`: the header file for `RPD\src\bbstep_PGM.cpp`.
- `RPD\src\cost.hpp`: the header file containing the link performance function that used to solve either lower-level traffic assignment problem or Problem (18).
- `RPD\src\path.hpp`: the header file containing the structure of the path.
- `RPD\src\graph.hpp`: the header file containing the structure of the graph.
- `RPD\src\io.hpp`: the header file to handle input/output streams.
- `RPD\src\utils.hpp`: the header file containing auxiliary functions.
- `RPD\src\path_equilibration.cpp`: the source codes for solving the lower-level traffic assignment problem with given $u$.
- `RPD\src\path_equilibration_rho.cpp`: the source codes for solving Problem (18) (i.e., step 1 of Algorithm 3) to update $v$.
- `RPD\src\bbstep_PGM.cpp`: the source codes for solving Problem (19) (i.e., step 2 of Algorithm 3) to update $u$.
- `RPD\main.cpp`: the source codes of Algorithm RPD.
