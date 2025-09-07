#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <math.h>

// Struct to represent a 3D coordinate point.
typedef struct {
    int x, y, z;
} Point;

// Converts a 1D index to 3D coordinates within a subdomain of size X × Y × Z.
Point pos_to_coords(int i, int X, int Y, int Z) {
    Point coords;
    coords.z = i / (X * Y);
    int remainder = i % (X * Y);
    coords.y = remainder / X;
    coords.x = remainder % X;
    return coords;
}

// Converts 3D coordinates to a 1D index within a subdomain of size X × Y × Z.
int coords_to_pos(Point coords, int X, int Y, int Z) {
    return coords.z * X * Y + coords.y * X + coords.x;
}

// Determines whether a point is a local maximum or minimum by comparing with local and halo neighbors.
void isLocalMaxMin(bool *isMin, bool *isMax, float *local_data, int col, int local_index, int NC, int NX, int NY, int NZ, int PX, int PY, int PZ,float *X_plus_ones_layer_recv,float *X_minus_ones_layer_recv,float *Y_plus_ones_layer_recv,float *Y_minus_ones_layer_recv,float *Z_plus_ones_layer_recv,float *Z_minus_ones_layer_recv,Point rankCoordinates){
    Point local_coordinates = pos_to_coords(local_index, NX/PX, NY/PY, NZ/PZ);
    
    *isMax = true;
    *isMin = true;
    
    if(((*isMin) || (*isMax)) && local_coordinates.x == 0 && rankCoordinates.x > 0){
        int neighbour_index = (NY/PY)*local_coordinates.z + local_coordinates.y;
        if(X_minus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(X_minus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }
    
    if(((*isMin) || (*isMax)) && local_coordinates.x == NX/PX-1 && rankCoordinates.x < PX-1){
        int neighbour_index = (NY/PY)*local_coordinates.z + local_coordinates.y;
        if(X_plus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(X_plus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }

    if(((*isMin) || (*isMax)) && local_coordinates.y == 0 && rankCoordinates.y > 0){
        int neighbour_index = (NX/PX)*local_coordinates.z + local_coordinates.x;
        if(Y_minus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(Y_minus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }
    
    if(((*isMin) || (*isMax)) && local_coordinates.y == NY/PY-1 && rankCoordinates.y < PY-1){
        int neighbour_index = (NX/PX)*local_coordinates.z + local_coordinates.x;
        if(Y_plus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(Y_plus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }
    
    if(((*isMin) || (*isMax)) && local_coordinates.z == 0 && rankCoordinates.z > 0){
        int neighbour_index = (NX/PX)*local_coordinates.y + local_coordinates.x;
        if(Z_minus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(Z_minus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }
    
    if(((*isMin) || (*isMax)) && local_coordinates.z == NZ/PZ-1 && rankCoordinates.z < PZ-1){
        int neighbour_index = (NX/PX)*local_coordinates.y + local_coordinates.x;
        if(Z_plus_ones_layer_recv[neighbour_index]>=local_data[local_index*NC + col])
        *isMax = false;
        if(Z_plus_ones_layer_recv[neighbour_index]<=local_data[local_index*NC + col])
        *isMin = false;
    }
    
    if((*isMax)
        && (local_coordinates.x == 0 || local_data[local_index*NC + col] > local_data[coords_to_pos((Point){local_coordinates.x-1,local_coordinates.y,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.x == NX/PX-1 || local_data[local_index*NC + col]> local_data[coords_to_pos((Point){local_coordinates.x+1,local_coordinates.y,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.y == 0 || local_data[local_index*NC + col]> local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y-1,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.y == NY/PY-1 || local_data[local_index*NC + col]> local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y+1,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.z == 0 || local_data[local_index*NC + col]> local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y,local_coordinates.z-1},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.z == NZ/PZ-1 || local_data[local_index*NC + col]> local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y,local_coordinates.z+1},NX/PX,NY/PY,NZ/PZ)*NC + col])
    ){
        *isMax = true;
    }
    else {
        *isMax = false;
    }

    if((*isMin)
        && (local_coordinates.x == 0 || local_data[local_index*NC + col] < local_data[coords_to_pos((Point){local_coordinates.x-1,local_coordinates.y,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.x == NX/PX-1 || local_data[local_index*NC + col]< local_data[coords_to_pos((Point){local_coordinates.x+1,local_coordinates.y,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.y == 0 || local_data[local_index*NC + col]< local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y-1,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.y == NY/PY-1 || local_data[local_index*NC + col]< local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y+1,local_coordinates.z},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.z == 0 || local_data[local_index*NC + col]< local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y,local_coordinates.z-1},NX/PX,NY/PY,NZ/PZ)*NC + col])
        && (local_coordinates.z == NZ/PZ-1 || local_data[local_index*NC + col] < local_data[coords_to_pos((Point){local_coordinates.x,local_coordinates.y,local_coordinates.z+1},NX/PX,NY/PY,NZ/PZ)*NC + col])
    ){
        *isMin = true;
    }
    else {
        *isMin = false;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 10) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <input_filename> <PX> <PY> <PZ> <NX> <NY> <NZ> <NC> <output_filename>\n", argv[0]);
            fprintf(stderr, "Error: Incorrect number of arguments\n");
        }
        MPI_Finalize();
        return -1;
    }
    
    char *input_filename = argv[1];
    int PX = atoi(argv[2]);
    int PY = atoi(argv[3]);
    int PZ = atoi(argv[4]);
    int NX = atoi(argv[5]);
    int NY = atoi(argv[6]);
    int NZ = atoi(argv[7]);
    int NC = atoi(argv[8]);
    char *output_filename = argv[9];
    
    
    if (PX * PY * PZ != size) {
        if (rank == 0) {
            fprintf(stderr, "Error: Number of processes (%d) does not match the product of PX, PY, and PZ (%d)\n", size, PX * PY * PZ);
            fprintf(stderr, "Error: PX*PY*PZ must equal the number of processes\n");
        }
        MPI_Finalize();
        return -1;
    }
    
    double time1, time2, time3, time4;
    
    Point rank_coordinates = pos_to_coords(rank, PX, PY, PZ);
    int elements_per_rank = (NX * NY * NZ) / (PX * PY * PZ);
    int NX_per_process = NX / PX;
    int NY_per_process = NY / PY;
    int NZ_per_process = NZ / PZ;
    
    float* local_data = (float*)malloc(NC * elements_per_rank * sizeof(float));

    MPI_Datatype file_subarray;
    int global_sizes[4] = {NZ, NY, NX, NC};  // File storage order (slowest to fastest: Z, Y, X, NC)
    int local_sizes[4] = {NZ_per_process, NY_per_process, NX_per_process, NC};
    int starts[4] = {rank_coordinates.z * NZ_per_process,  rank_coordinates.y* NY_per_process, rank_coordinates.x * NX_per_process, 0};

    MPI_Type_create_subarray(
        4,                      // 4D array (Z, Y, X, NC)
        global_sizes,           // Global size
        local_sizes,            // Local size
        starts,                 // Starting indices
        MPI_ORDER_C,            // C-style row-major (X fastest, Z slowest)
        MPI_FLOAT,              // Data type
        &file_subarray
    );
    MPI_Type_commit(&file_subarray);

    time1 = MPI_Wtime();

    MPI_File fh;
    MPI_File_open(
        MPI_COMM_WORLD,
        input_filename,
        MPI_MODE_RDONLY,
        MPI_INFO_NULL,
        &fh
    );

    MPI_File_set_view(
        fh,
        0,                      // Offset (bytes)
        MPI_FLOAT,              // Displacement type
        file_subarray,          // File view (4D subarray)
        "native",               // Data representation
        MPI_INFO_NULL
    );

    MPI_File_read_all(
        fh,
        local_data,
        elements_per_rank * NC,
        MPI_FLOAT,
        MPI_STATUS_IGNORE
    );


    time2 = MPI_Wtime();


    float** X_plus_ones_layer_send = (float**)malloc(NC * sizeof(float*));
    float** X_minus_ones_layer_send = (float**)malloc(NC * sizeof(float*));
    float** Y_plus_ones_layer_send = (float**)malloc(NC * sizeof(float*));
    float** Y_minus_ones_layer_send = (float**)malloc(NC * sizeof(float*));
    float** Z_plus_ones_layer_send = (float**)malloc(NC * sizeof(float*));
    float** Z_minus_ones_layer_send = (float**)malloc(NC * sizeof(float*));

    float** X_plus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));
    float** X_minus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));
    float** Y_plus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));
    float** Y_minus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));
    float** Z_plus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));
    float** Z_minus_ones_layer_recv = (float**)malloc(NC * sizeof(float*));

    // Allocate memory for the layers
    for(int i = 0; i < NC; i++) {
        X_plus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NX_per_process * sizeof(float));
        X_minus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NX_per_process * sizeof(float));
        Y_plus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NY_per_process * sizeof(float));
        Y_minus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NY_per_process * sizeof(float));
        Z_plus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NZ_per_process * sizeof(float));
        Z_minus_ones_layer_send[i] = (float*)malloc(elements_per_rank/NZ_per_process * sizeof(float));
        
        X_plus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NX_per_process * sizeof(float));
        X_minus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NX_per_process * sizeof(float));
        Y_plus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NY_per_process * sizeof(float));
        Y_minus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NY_per_process * sizeof(float));
        Z_plus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NZ_per_process * sizeof(float));
        Z_minus_ones_layer_recv[i] = (float*)malloc(elements_per_rank/NZ_per_process * sizeof(float));
    }

        
    MPI_Request *r1 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r2 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r3 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r4 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r5 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r6 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r7 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r8 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r9 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r10 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r11 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));
    MPI_Request *r12 = (MPI_Request*)malloc(NC * sizeof(MPI_Request));

    for(int i=0;i<NC;i++){
        r1[i] = MPI_REQUEST_NULL;
        r2[i] = MPI_REQUEST_NULL;
        r3[i] = MPI_REQUEST_NULL;
        r4[i] = MPI_REQUEST_NULL;
        r5[i] = MPI_REQUEST_NULL;
        r6[i] = MPI_REQUEST_NULL;
        r7[i] = MPI_REQUEST_NULL;
        r8[i] = MPI_REQUEST_NULL;
        r9[i] = MPI_REQUEST_NULL;
        r10[i] = MPI_REQUEST_NULL;
        r11[i] = MPI_REQUEST_NULL;
        r12[i] = MPI_REQUEST_NULL;
    }
    
    // communication with X-1 rank
    if(rank_coordinates.x > 0) {
        Point X_minus_one_coords = {rank_coordinates.x - 1, rank_coordinates.y, rank_coordinates.z};
        int X_minus_one_rank = coords_to_pos(X_minus_one_coords, PX, PY, PZ);
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int i = 0; i < elements_per_rank; i+= NX_per_process){
                X_minus_ones_layer_send[t][index++] = local_data[i*NC + t];
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(X_minus_ones_layer_send[t], elements_per_rank/NX_per_process, MPI_FLOAT, X_minus_one_rank, rank*NC + t, MPI_COMM_WORLD, &r1[t]);
        }

        for(int t =0;t<NC;t++){
            MPI_Irecv(X_minus_ones_layer_recv[t], elements_per_rank/NX_per_process, MPI_FLOAT, X_minus_one_rank, X_minus_one_rank*NC + t + size*NC, MPI_COMM_WORLD, &r2[t]);
        }
        
        // printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received X-1 layer from rank %d for testcase %d:\n", rank, X_minus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NX_per_process; j++) {
        //         printf("%f ", X_minus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }

    // communication with X+1 rank
    if(rank_coordinates.x < PX - 1) {
        
        Point X_plus_one_coords = {rank_coordinates.x + 1, rank_coordinates.y, rank_coordinates.z};
        int X_plus_one_rank = coords_to_pos(X_plus_one_coords, PX, PY, PZ);
        
        for(int t= 0;t<NC;t++){
            MPI_Irecv(X_plus_ones_layer_recv[t], elements_per_rank/NX_per_process, MPI_FLOAT, X_plus_one_rank, X_plus_one_rank*NC + t, MPI_COMM_WORLD, &r3[t]);
        }
        
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int i = NX_per_process-1; i < elements_per_rank; i+= NX_per_process){
                X_plus_ones_layer_send[t][index++] = local_data[i*NC + t];
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(X_plus_ones_layer_send[t], elements_per_rank/NX_per_process, MPI_FLOAT, X_plus_one_rank, rank*NC + t + size*NC, MPI_COMM_WORLD, &r4[t]);
        }


        //printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received X+1 layer from rank %d for testcase %d:\n", rank, X_plus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NX_per_process; j++) {
        //         printf("%f ", X_plus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }

    // communication with Y-1 rank
    if(rank_coordinates.y > 0) {
        Point Y_minus_one_coords = {rank_coordinates.x, rank_coordinates.y - 1, rank_coordinates.z};
        int Y_minus_one_rank = coords_to_pos(Y_minus_one_coords, PX, PY, PZ);
        
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int z = 0;z<NZ_per_process;z++){
                for(int x=0;x<NX_per_process;x++){
                    Y_minus_ones_layer_send[t][index++] = local_data[(z*NX_per_process*NY_per_process + x)*NC + t];
                }
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(Y_minus_ones_layer_send[t], elements_per_rank/NY_per_process, MPI_FLOAT, Y_minus_one_rank, rank*NC + t + size*NC*2, MPI_COMM_WORLD, &r5[t]);
        }

        for(int t =0;t<NC;t++){
            MPI_Irecv(Y_minus_ones_layer_recv[t], elements_per_rank/NY_per_process, MPI_FLOAT, Y_minus_one_rank, Y_minus_one_rank*NC + t + size*NC*3, MPI_COMM_WORLD, &r6[t]);
        }
        
        // printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received Y-1 layer from rank %d for testcase %d:\n", rank, Y_minus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NY_per_process; j++) {
        //         printf("%f ", Y_minus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }


    // communication with Y+1 rank
    if(rank_coordinates.y < PY - 1) {
        Point Y_plus_one_coords = {rank_coordinates.x, rank_coordinates.y + 1, rank_coordinates.z};
        int Y_plus_one_rank = coords_to_pos(Y_plus_one_coords, PX, PY, PZ);
        
        for(int t= 0;t<NC;t++){
            MPI_Irecv(Y_plus_ones_layer_recv[t], elements_per_rank/NY_per_process, MPI_FLOAT, Y_plus_one_rank, Y_plus_one_rank*NC + t + size*NC*2, MPI_COMM_WORLD, &r7[t]);
        }
        
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int z = 0;z<NZ_per_process;z++){
                for(int x = 0;x<NX_per_process;x++){
                    Y_plus_ones_layer_send[t][index++] = local_data[(z*NX_per_process*NY_per_process + (NY_per_process-1)*NX_per_process + x)*NC + t];
                }
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(Y_plus_ones_layer_send[t], elements_per_rank/NY_per_process, MPI_FLOAT, Y_plus_one_rank, rank*NC + t + size*NC*3, MPI_COMM_WORLD, &r8[t]);
        }
        //printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received Y+1 layer from rank %d for testcase %d:\n", rank, Y_plus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NY_per_process; j++) {
        //         printf("%f ", Y_plus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }


    // communication with Z-1 rank
    if(rank_coordinates.z > 0) {
        Point Z_minus_one_coords = {rank_coordinates.x, rank_coordinates.y, rank_coordinates.z - 1};
        int Z_minus_one_rank = coords_to_pos(Z_minus_one_coords, PX, PY, PZ);
        
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int i = 0; i<NX_per_process*NY_per_process; i++){
                Z_minus_ones_layer_send[t][index++] = local_data[i*NC + t];
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(Z_minus_ones_layer_send[t], elements_per_rank/NZ_per_process, MPI_FLOAT, Z_minus_one_rank, rank*NC + t + size*NC*4, MPI_COMM_WORLD, &r9[t]);
        }

        for(int t =0;t<NC;t++){
            MPI_Irecv(Z_minus_ones_layer_recv[t], elements_per_rank/NZ_per_process, MPI_FLOAT, Z_minus_one_rank, Z_minus_one_rank*NC + t + size*NC*5, MPI_COMM_WORLD, &r10[t]);
        }
        
        // printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received Z-1 layer from rank %d for testcase %d:\n", rank, Z_minus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NZ_per_process; j++) {
        //         printf("%f ", Z_minus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }

    // communication with Z+1 rank
    if(rank_coordinates.z < PZ - 1) {
        Point Z_plus_one_coords = {rank_coordinates.x, rank_coordinates.y, rank_coordinates.z + 1};
        int Z_plus_one_rank = coords_to_pos(Z_plus_one_coords, PX, PY, PZ);
        
        for(int t= 0;t<NC;t++){
            MPI_Irecv(Z_plus_ones_layer_recv[t], elements_per_rank/NZ_per_process, MPI_FLOAT, Z_plus_one_rank, Z_plus_one_rank*NC + t + size*NC*4, MPI_COMM_WORLD, &r11[t]);
        }
        
        for(int t= 0;t<NC;t++){
            int index = 0;
            for(int i = 0; i < NX_per_process*NY_per_process; i++){
                Z_plus_ones_layer_send[t][index++] = local_data[(i+(NZ_per_process-1)*NX_per_process*NY_per_process)*NC + t];
            }
        }
        
        for(int t =0;t<NC;t++){
            MPI_Isend(Z_plus_ones_layer_send[t], elements_per_rank/NZ_per_process, MPI_FLOAT, Z_plus_one_rank, rank*NC + t + size*NC*5, MPI_COMM_WORLD, &r12[t]);
        }
        
        //printing the received layer
        // for(int i = 0; i < NC; i++) {
        //     printf("Rank %d received Z+1 layer from rank %d for testcase %d:\n", rank, Z_plus_one_rank, i);
        //     for(int j = 0; j < elements_per_rank/NZ_per_process; j++) {
        //         printf("%f ", Z_plus_ones_layer_recv[i][j]);
        //     }
        //     printf("\n");
        // }
    }
    
    // Checking for the local minimum and maximums
    int local_minimaCount[NC], local_maximaCount[NC];
    float local_maximum[NC], local_minimum[NC];
    for(int i = 0; i < NC; i++) {
        local_maximum[i] = local_data[i];
        local_minimum[i] = local_data[i];
        local_minimaCount[i] = 0;
        local_maximaCount[i] = 0;
    }

    MPI_Waitall(NC, r1, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r2, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r3, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r4, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r5, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r6, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r7, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r8, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r9, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r10, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r11, MPI_STATUSES_IGNORE);
    MPI_Waitall(NC, r12, MPI_STATUSES_IGNORE);

    for(int t=0;t<NC;t++){
        for(int i=0;i<elements_per_rank;i++){
            bool isMin = false, isMax = false;
            isLocalMaxMin(&isMin, &isMax, local_data, t, i, NC, NX, NY, NZ, PX, PY, PZ, X_plus_ones_layer_recv[t], X_minus_ones_layer_recv[t], Y_plus_ones_layer_recv[t], Y_minus_ones_layer_recv[t], Z_plus_ones_layer_recv[t], Z_minus_ones_layer_recv[t], rank_coordinates);
            if(isMin){
                local_minimaCount[t]++;
            }
            if(isMax){
                local_maximaCount[t]++;
            }
            
            if(local_data[i*NC + t] > local_maximum[t]){
                local_maximum[t] = local_data[i*NC + t];
            }
            if(local_data[i*NC + t] < local_minimum[t]){
                local_minimum[t] = local_data[i*NC + t];
            }
        }
    }
    
    int global_minimaCount[NC], global_maximaCount[NC];
    float global_maximum[NC], global_minimum[NC];
    MPI_Reduce(local_minimaCount, global_minimaCount, NC, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_maximaCount, global_maximaCount, NC, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_maximum, global_maximum, NC, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_minimum, global_minimum, NC, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    
    time3 = MPI_Wtime();



    double read_time = time2 - time1;
    double main_code_time = time3 - time2;
    double total_time = time3 - time1;
    double read_time_max, main_code_time_max, total_time_max;
    MPI_Reduce(&read_time, &read_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&main_code_time, &main_code_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &total_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if(!rank){
        FILE *output_file = fopen(output_filename, "w");
        if (output_file == NULL) {
            printf("Error opening output file\n");
            MPI_Finalize();
            return -1;
        }

        for(int i = 0; i < NC; i++) {
            fprintf(output_file, "(%d, %d), ", global_minimaCount[i], global_maximaCount[i]);
        }
        fprintf(output_file, "\n");

        for(int i = 0; i < NC; i++) {
            fprintf(output_file, "(%.4f, %.4f), ", global_minimum[i], global_maximum[i]);
        }
        fprintf(output_file, "\n");

        fprintf(output_file, "%lf, %lf, %lf\n", read_time_max, main_code_time_max, total_time_max);

        fclose(output_file);
        printf("Output written to %s\n", output_filename);
    }


    free(local_data);

    MPI_Finalize();
    return 0;
}