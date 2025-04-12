#include "bits/stdc++.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <mpi.h>

using namespace std;

typedef long long ll;

int PX, PY, PZ, nx, ny, nz, t, processes;
double data_dist_time, compute_time, total_time;
vector<int> minima, maxima, total_minima, total_maxima;

// utilities
int get_data(int nx, int ny, int nz, int t, string &filename, vector<float> &data)
{
    ifstream infile(filename);
    if (!infile.is_open())
    {
        cerr << "Could not open file." << endl;
        return 1;
    }

    string line;
    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                if (!getline(infile, line))
                {
                    cerr << "Unexpected end of file." << endl;
                    return 1;
                }
                stringstream ss(line);
                for (int i = 0; i < t; ++i)
                {
                    ll idx = x * (ny * nz * t) + y * (nz * t) + z * t + i;
                    if (!(ss >> data[idx]))
                    {
                        cerr << "Error reading data." << endl;
                        return 1;
                    }
                }
            }
        }
    }
    infile.close();
    return 0;
}

int get_data_binary(int nx, int ny, int nz, int t, string &filename, vector<float> &data)
{
    cout << "Trying to open file: " << filename << endl;
    ifstream infile(filename, ios::binary);
    if (!infile.is_open())
    {
        cerr << "Could not open binary file." << endl;
        return 1;
    }

    size_t total_elements = nx * ny * nz * t;
    vector<float> raw_data(total_elements);
    // data.resize(total_elements);

    // Read all data at once
    infile.read(reinterpret_cast<char *>(raw_data.data()), total_elements * sizeof(float));
    if (!infile)
    {
        cerr << "Error reading binary data." << endl;
        return 1;
    }

    // Rearranging according to your custom layout
    size_t linear_idx = 0;
    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                for (int i = 0; i < t; ++i)
                {
                    size_t idx = x * (ny * nz * t) + y * (nz * t) + z * t + i;
                    data[idx] = raw_data[linear_idx++];
                }
            }
        }
    }

    infile.close();
    return 0;
}

tuple<int, int, int, int, int, int> get_bounding_box(int rank)
{

    int x_partition_len, y_partition_len, z_partition_len, x_req = 0, y_req = 0, z_req = 0;

    x_partition_len = nx / PX;
    x_req = nx % PX;
    y_partition_len = ny / PY;
    y_req = ny % PY;
    z_partition_len = nz / PZ;
    z_req = nz % PZ;

    int x_pos = rank % PX, y_pos = (rank / PX) % PY, z_pos = rank / (PX * PY);

    int start_X = x_pos * (x_partition_len + 1), end_X = min(start_X + x_partition_len + 1, nx);
    int start_Y = y_pos * (y_partition_len + 1), end_Y = min(start_Y + y_partition_len + 1, ny);
    int start_Z = z_pos * (z_partition_len + 1), end_Z = min(start_Z + z_partition_len + 1, nz);

    if (x_pos >= x_req)
    {
        start_X = (x_req * (x_partition_len + 1)) + ((x_pos - x_req) * x_partition_len);
        end_X = min(start_X + x_partition_len, nx);
    }

    if (y_pos >= y_req)
    {
        start_Y = (y_req * (y_partition_len + 1)) + ((y_pos - y_req) * y_partition_len);
        end_Y = min(start_Y + y_partition_len, ny);
    }
    if (z_pos >= z_req)
    {
        start_Z = (z_req * (z_partition_len + 1)) + ((z_pos - z_req) * z_partition_len);
        end_Z = min(start_Z + z_partition_len, nz);
    }
    return {start_X, end_X, start_Y, end_Y, start_Z, end_Z};
}

tuple<int, int> minimax(int x, int y, int z, int time_idx, vector<vector<vector<vector<float>>>> &grid)
{
    int a[] = {1, -1, 0, 0, 0, 0};
    int b[] = {0, 0, 1, -1, 0, 0};
    int c[] = {0, 0, 0, 0, 1, -1};

    int is_min = true, is_max = true;
    float curr = grid[x][y][z][time_idx];

    for (int i = 0; i < 6; i++)
    {
        int new_x = x + a[i];
        int new_y = y + b[i];
        int new_z = z + c[i];

        if (new_x < 0 || new_y < 0 || new_z < 0 || new_x == grid.size() || new_y == grid[0].size() || new_z == grid[0][0].size())
            continue;

        if (grid[new_x][new_y][new_z][time_idx] <= curr)
            is_min = false;

        if (grid[new_x][new_y][new_z][time_idx] >= curr)
            is_max = false;
    }

    return {is_min, is_max};
}

void compute(vector<vector<vector<vector<float>>>> &grid, int rank)
{
    auto [start_X, end_X, start_Y, end_Y, start_Z, end_Z] = get_bounding_box(rank);

    int ex = end_X - start_X, ey = end_Y - start_Y, ez = end_Z - start_Z;

    ex += (start_X != 0);
    ey += (start_Y != 0);
    ez += (start_Z != 0);

    int sx = start_X != 0, sy = start_Y != 0, sz = start_Z != 0;
    // int local_minima = 0, local_maxima = 0;
    for (int time_idx = 0; time_idx < t; time_idx++)
    {
        int local_minima = 0, local_maxima = 0;
        for (int x = sx; x < ex; x++)
        {
            for (int y = sy; y < ey; y++)
            {
                for (int z = sz; z < ez; z++)
                {
                    auto [is_min, is_max] = minimax(x, y, z, time_idx, grid);
                    local_minima += is_min;
                    local_maxima += is_max;
                }
            }
        }
        minima[time_idx] = local_minima;
        maxima[time_idx] = local_maxima;
    }
    // return {local_minima, local_maxima};
    return;
}

int main(int argc, char *argv[])
{
    if (argc != 10)
    {
        cerr << "Usage: " << argv[0] << " <Dataset> <PX> <PY> <PZ> <NX> <NY> <NZ> <NC> <Output File>" << endl;
        exit(1);
    }

    string filename = argv[1];
    PX = stoi(argv[2]);
    PY = stoi(argv[3]);
    PZ = stoi(argv[4]);
    nx = stoi(argv[5]); // number of grid points along X
    ny = stoi(argv[6]); // number of grid points along Y
    nz = stoi(argv[7]); // number of grid points along Z
    t = stoi(argv[8]);  // number of time stamps
    string output_filename = argv[9];

    processes = PX * PY * PZ;

    MPI_Init(&argc, &argv);
    double t_1 = MPI_Wtime(), compute_start, compute_end, total_compute = 0, comm_start, comm_end, ray_tracing_time;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    vector<float> data;
    vector<vector<vector<vector<float>>>> grid;
    vector<float> global_minima_val(t, INT_MAX), global_maxima_val(t, INT_MIN);
    if (rank == 0)
    {
        // vector<float> data;
        cout << nx << " " << ny << " " << nz << " " << t << endl;
        ll total_size = (ll)nx * ny * nz * t;
        try
        {
            data.resize(total_size);
        }
        catch (const std::bad_alloc &e)
        {
            cerr << "Memory allocation failed: " << e.what() << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (get_data_binary(nx, ny, nz, t, filename, data))
        {
            cerr << "failed loading" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < 20; i++)
            cout << data[i] << endl;

        for (int i = 1; i < processes; i++)
        {
            auto [start_X, end_X, start_Y, end_Y, start_Z, end_Z] = get_bounding_box(i);
            int x_siz = end_X - start_X, y_siz = end_Y - start_Y, z_siz = end_Z - start_Z;
            x_siz += (start_X != 0) + (end_X != nx);
            y_siz += (start_Y != 0) + (end_Y != ny);
            z_siz += (start_Z != 0) + (end_Z != nz);

            int siz = x_siz * y_siz * z_siz * t, temp_idx = 0;
            vector<float> temp(siz);

            for (int time_idx = 0; time_idx < t; time_idx++)
            {
                for (int z = max(0, start_Z - 1); z < min(nz, end_Z + 1); z++)
                {
                    for (int y = max(0, start_Y - 1); y < min(ny, end_Y + 1); y++)
                    {
                        for (int x = max(0, start_X - 1); x < min(nx, end_X + 1); x++)
                        {
                            ll data_idx = x * (ny * nz * t) + y * (nz * t) + z * t + time_idx;
                            // cout << data[data_idx] << " ";
                            temp[temp_idx++] = data[data_idx];
                            global_minima_val[time_idx] = min(global_minima_val[time_idx], data[data_idx]);
                            global_maxima_val[time_idx] = max(global_maxima_val[time_idx], data[data_idx]);
                        }
                    }
                }
            }

            MPI_Send(&temp[0], siz, MPI_FLOAT, i, 99, MPI_COMM_WORLD);
            temp.clear();
        }

        auto [start_X, end_X, start_Y, end_Y, start_Z, end_Z] = get_bounding_box(rank);
        int x_siz = end_X - start_X, y_siz = end_Y - start_Y, z_siz = end_Z - start_Z;
        x_siz += (start_X != 0) + (end_X != nx);
        y_siz += (start_Y != 0) + (end_Y != ny);
        z_siz += (start_Z != 0) + (end_Z != nz);
        grid = vector<vector<vector<vector<float>>>>(x_siz, vector<vector<vector<float>>>(y_siz, vector<vector<float>>(z_siz, vector<float>(t))));

        for (int x = max(0, start_X - 1); x < min(nx, end_X + 1); x++)
        {
            for (int y = max(0, start_Y - 1); y < min(ny, end_Y + 1); y++)
            {
                for (int z = max(0, start_Z - 1); z < min(nz, end_Z + 1); z++)
                {
                    for (int time_idx = 0; time_idx < t; ++time_idx)
                    {
                        ll data_idx = x * (ny * nz * t) + y * (nz * t) + z * t + time_idx;
                        grid[x - max(0, start_X - 1)][y - max(0, start_Y - 1)][z - max(0, start_Z - 1)][time_idx] = data[data_idx];
                    }
                }
            }
        }
        data.clear();
    }
    else
    {
        auto [start_X, end_X, start_Y, end_Y, start_Z, end_Z] = get_bounding_box(rank);
        int x_siz = end_X - start_X, y_siz = end_Y - start_Y, z_siz = end_Z - start_Z;
        x_siz += (start_X != 0) + (end_X != nx);
        y_siz += (start_Y != 0) + (end_Y != ny);
        z_siz += (start_Z != 0) + (end_Z != nz);

        int siz = x_siz * y_siz * z_siz * t, temp_idx = 0;
        data.resize(siz);

        MPI_Recv(&data[0], siz, MPI_FLOAT, 0, 99, MPI_COMM_WORLD, &status);

        grid = vector<vector<vector<vector<float>>>>(x_siz, vector<vector<vector<float>>>(y_siz, vector<vector<float>>(z_siz, vector<float>(t))));
        int data_idx = 0;
        for (int time_idx = 0; time_idx < t; ++time_idx)
        {
            for (int z = 0; z < z_siz; z++)
            {
                for (int y = 0; y < y_siz; y++)
                {
                    for (int x = 0; x < x_siz; x++)
                    {
                        grid[x][y][z][time_idx] = data[data_idx++];
                    }
                }
            }
        }

        data.clear();
    }
    double t_2 = MPI_Wtime();

    minima.resize(t, 0);
    maxima.resize(t, 0);
    compute(grid, rank);
    double t_3 = MPI_Wtime();
    double local_data_dist_time = t_2 - t_1, local_total_time = t_3 - t_1, local_compute_time = t_3 - t_2;
    if (rank == 0)
    {
        total_minima.resize(t), total_maxima.resize(t);
    }

    MPI_Reduce(minima.data(), total_minima.data(), t, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(maxima.data(), total_maxima.data(), t, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_data_dist_time, &data_dist_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_total_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_compute_time, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        string out = "output_";
        out += (to_string(nx) + "_");
        out += (to_string(ny) + "_");
        out += (to_string(nz) + "_");
        out += to_string(t);
        out += ".txt";
        ofstream outfile(out);

        for (int i = 0; i < t; i++)
        {
            outfile << "( " << total_minima[i] << " , " << total_maxima[i] << " )" << endl;
        }
        outfile.precision(6);
        outfile << std::fixed;
        for (int i = 0; i < t; i++)
        {
            outfile << "( " << global_minima_val[i] << " , " << global_maxima_val[i] << " )" << endl;
        }

        outfile << "compute time -> " << compute_time << endl;
        outfile << "data read and dist time -> " << data_dist_time << endl;
        outfile << "total time -> " << total_time << endl;

        outfile.close();
    }
    MPI_Finalize();
}
