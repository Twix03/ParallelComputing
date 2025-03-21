# ParallelComputing
This repo will contain the assignment solutions for CS633.

### Running the code 

0. Make sure to install mpi and check for required libraries for c++ ie {iostream, fstream, vector, cmath, sstream, algorithm, <filesystem>, <mpi.h>}.

1. To run on _x_ no. of cores , on the terminal execute :-
   
   ```
   mpicxx -std=c++20 -o executable  solution.cpp
   ```
   To run on single host :->
   ```
   mpirun -np 16 ./executable data_64_64_64_7.txt 4 4 1 64 64 64 7 output_64_64_64_3.txt
   ```
   To run on multiple hosts :->  
   ```
   mpirun -np 8 -f hostfile ./executable data_64_64_64_3.txt 2 2 2 64 64 64 3 output_64_64_64_3.txt
   ```

# Input (Nine Arguments)

1. **Dataset**: Input file name (e.g., `data.txt`)

2. **PX**: Number of processes in the X-dimension

3. **PY**: Number of processes in the Y-dimension

4. **PZ**: Number of processes in the Z-dimension

5. **NX**: Number of grid points in the X-dimension

6. **NY**: Number of grid points in the Y-dimension

7. **NZ**: Number of grid points in the Z-dimension

8. **NC**: Number of columns (number of time steps)

9. **Output File Name**: Preferably in the format `output_NX_NY_NZ_NC.txt`  
   - **Note**: Write the output to this file from rank 0.  

---

### Important Notes
- The number of processes `P` is given by:  
  \[
  P = PX \times PY \times PZ
  \]
- The number of processes (`P`) specified during execution should match this value.  
  - Use `-n` or `--np` (depending on the system) to specify the number of processes while running the application.
