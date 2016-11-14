#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <mpi.h>

/***********************************************************
  Helper functions
***********************************************************/

//For exiting on error condition
void die(int lineNo);

//For trackinng execution
long long wallClockTime();


/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue );

void freeSquareMatrix( char** );

void printSquareMatrix( char**, int size );


/***********************************************************
   World  related functions
***********************************************************/

#define ALIVE 'X'
#define DEAD 'O'

char** readWorldFromFile( char* fname, int* size );

int countNeighbours(char** world, int row, int col);

void evolveWorld(char** curWorld, char** nextWorld, int size, int iter, MPI_Comm cartesian_comm);


/***********************************************************
   Simple circular linked list for match records
***********************************************************/
typedef struct BMSTRUCT {
    int iteration, row, col, rotation;
} BASEMATCH;

int compareBaseMatch(const void* a, const void* b);

typedef struct MSTRUCT {
    int iteration, row, col, rotation;
    struct MSTRUCT *next;
} MATCH;

typedef struct {
    int nItem;
    MATCH* tail;
} MATCHLIST;

MATCHLIST* newList();

void deleteList(MATCHLIST*);

void insertEnd(MATCHLIST*, int, int, int, int);

void printList(MATCHLIST*);


/***********************************************************
   Search related functions
***********************************************************/

//Using the compass direction to indicate the rotation of pattern
#define N 0 //no rotation
#define E 1 //90 degree clockwise
#define S 2 //180 degree clockwise
#define W 3 //90 degree anti-clockwise

char** readPatternFromFile( char* fname, int* size );

void rotate90(char** current, char** rotated, int size);

void searchPatterns(char** world, int wSize, int iteration,
        char** patterns[4], int pSize, MATCHLIST* list);

void searchSinglePattern(char** world, int wSize, int interation,
        char** pattern, int pSize, int rotation, MATCHLIST* list);


/***********************************************************
   Cartesian block-wise related definitions
***********************************************************/

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
#define UPLEFT 4
#define UPRIGHT 5
#define DOWNLEFT 6
#define DOWNRIGHT 7


/***********************************************************
   Parallel related variables and functions
***********************************************************/

int rank, numtasks;
int evolverProcessors, finderProcessors;
int rowCount, colCount;

int rowSize, colSize;
int rowStart, rowEnd, colStart, colEnd;
int currRowSize, currColSize;

int cartRank, cartNumtasks;
int neighbours[8], coordinate[2];

int getCorrespondingEvolverOrFinderRank(int rank);
int getTag(int iter, int senderRank);
int getDiagonalCartesianNeighbour(int dir);
void getWorldFromEvolver(char** curWorld, int size, int iter, int patternSize);
void sendWorldToFinder(char ** curW, int size, int iter);


/***********************************************************
   Main function
***********************************************************/

int main( int argc, char** argv)
{
    char **curW, **nextW, **temp;
    char **patterns[4];
    int worldSize, patternSize;
    int dir, iterations, iter;
    long long before, after;
    MATCHLIST*list;

    if (argc < 4 ){
        fprintf(stderr,
            "Usage: %s <world file> <Iterations> <pattern file>\n", argv[0]);
        exit(1);
    }

    curW = readWorldFromFile(argv[1], &worldSize);
    nextW = allocateSquareMatrix(worldSize+2, DEAD);
    iterations = atoi(argv[2]);

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    if (numtasks < 4 || numtasks % 4 != 0) {
        fprintf(stderr,
            "Number of processors should be >= 4 and divisible by 4.\n");
        exit(0);
    }
    evolverProcessors = numtasks / 2;
    finderProcessors = numtasks / 2;
    rowCount = (int)ceil(sqrt(evolverProcessors));
    colCount = evolverProcessors / rowCount;

    // Create 2 communication groups -- 1 for evolving one generation to the
    // next, 1 for pattern finding.
    MPI_Group original_group, partition_group;

    int evolverRanks[evolverProcessors], finderRanks[finderProcessors];
    for (iter = 0; iter < evolverProcessors; iter++) {
        evolverRanks[iter] = iter;
    }
    for (iter = 0; iter < finderProcessors; iter++) {
        finderRanks[iter] = iter + evolverProcessors;
    }

    MPI_Comm_group(MPI_COMM_WORLD, &original_group);
    if (rank < evolverProcessors) {
        MPI_Group_incl(original_group, evolverProcessors, evolverRanks, &partition_group);
    }
    else {
        MPI_Group_incl(original_group, evolverProcessors, finderRanks, &partition_group);
    }

    // Create communication for partition within world communication.
    MPI_Comm partition_communication;
    MPI_Comm_create(MPI_COMM_WORLD, partition_group, &partition_communication);

    // Create cartesian within partition communication.
    int nDimension = 2;
    int dimension[2] = {rowCount, colCount};
    int periods[2] = {0,0};
    MPI_Comm cartesian_communication;
    MPI_Cart_create(partition_communication, nDimension, dimension, periods, 0,
        &cartesian_communication);
    MPI_Comm_rank(cartesian_communication, &cartRank);
    MPI_Cart_coords(cartesian_communication, cartRank, nDimension, coordinate);
    MPI_Cart_shift(cartesian_communication, 0, 1,
        &neighbours[UP], &neighbours[DOWN]);
    MPI_Cart_shift(cartesian_communication, 1, 1,
        &neighbours[LEFT], &neighbours[RIGHT]);

    // MPI Cart does not support diagonal neighbours.
    for (iter = UPLEFT; iter <= DOWNRIGHT; iter++) {
        neighbours[iter] = getDiagonalCartesianNeighbour(iter);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Divide cells into processes.
    rowSize = worldSize / rowCount;
    colSize = worldSize / colCount;
    rowStart = 1 + (rowSize * coordinate[0]);
    rowEnd = coordinate[0] == (rowCount - 1)
                ? worldSize
                : rowStart + rowSize - 1;
    colStart = 1 + (colSize * coordinate[1]);
    colEnd = coordinate[1] == (colCount - 1)
                ? worldSize
                : colStart + colSize - 1;
    currRowSize = rowEnd - rowStart + 1;
    currColSize = colEnd - colStart + 1;

    patterns[N] = readPatternFromFile(argv[3], &patternSize);
    for (dir = E; dir <= W; dir++){
        patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
        rotate90(patterns[dir-1], patterns[dir], patternSize);
    }

    if (currRowSize < patternSize || currColSize < patternSize) {
        printf("World is too small while processor is too many. Increase size of world or reduce processor number.\n");
        exit(0);
    }

    if (rank == 0) {
        printf("World Size = %d\n", worldSize);
        printf("Iterations = %d\n", iterations);
        printf("Pattern size = %d\n", patternSize);
    }

#ifdef DEBUG
    printSquareMatrix(patterns[N], patternSize);
    printSquareMatrix(patterns[E], patternSize);
    printSquareMatrix(patterns[S], patternSize);
    printSquareMatrix(patterns[W], patternSize);
#endif

    //Start timer
    before = wallClockTime();

    //Actual work start
    list = newList();

    for (iter = 0; iter < iterations; iter++){
        if (rank < evolverProcessors) {
            evolveWorld( curW, nextW, worldSize, iter, cartesian_communication);
            temp = curW;
            curW = nextW;
            freeSquareMatrix(temp);
            nextW = allocateSquareMatrix(worldSize + 2, DEAD);

            sendWorldToFinder(curW, worldSize, iter);
        }
        else {
            searchPatterns(curW, worldSize, iter, patterns, patternSize, list);

            getWorldFromEvolver(curW, worldSize, iter, patternSize);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    BASEMATCH* listArray = (BASEMATCH*)malloc(sizeof(BASEMATCH) * list->nItem);
    MATCH* current = list->tail;
    for (iter = 0; iter < list->nItem; iter++) {
        BASEMATCH temp;
        temp.iteration = current->iteration;
        temp.row = current->row;
        temp.col = current->col;
        temp.rotation = current->rotation;
        listArray[iter] = temp;
        current = current->next;
    }

    // Create MPI data type
    const int basic_match_items = 4;
    int basic_match_blocklengths[4] = {1,1,1,1};
    MPI_Datatype basic_match_types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype mpi_basic_match_type;
    MPI_Aint mpi_basic_match_offsets[4];
    mpi_basic_match_offsets[0] = offsetof(BASEMATCH, iteration);
    mpi_basic_match_offsets[1] = offsetof(BASEMATCH, row);
    mpi_basic_match_offsets[2] = offsetof(BASEMATCH, col);
    mpi_basic_match_offsets[3] = offsetof(BASEMATCH, rotation);
    MPI_Type_create_struct(basic_match_items, basic_match_blocklengths, mpi_basic_match_offsets, basic_match_types, &mpi_basic_match_type);
    MPI_Type_commit(&mpi_basic_match_type);

    int* sizeList = (int*)malloc(sizeof(int) * numtasks);
    MPI_Gather(&(list->nItem), 1, MPI_INT, sizeList, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int* offsetList = (int*)malloc(sizeof(int) * numtasks);
    int startOffset = 0;
    if (rank == 0) {
        for (iter = 0; iter < numtasks; iter++) {
            offsetList[iter] = startOffset;
            startOffset += sizeList[iter];
        }
    }

    BASEMATCH* combinedList = (BASEMATCH*)malloc(startOffset * sizeof(BASEMATCH));
    MPI_Gatherv(listArray, list->nItem, mpi_basic_match_type, combinedList, sizeList, offsetList, mpi_basic_match_type, 0, MPI_COMM_WORLD );
    qsort(combinedList, startOffset, sizeof(BASEMATCH), compareBaseMatch);
    if (rank == 0) {
        printf("List size = %d\n", startOffset);
        for (iter = 0; iter < startOffset; iter++) {
            printf("%d:%d:%d:%d\n",
                    combinedList[iter].iteration,
                    combinedList[iter].row,
                    combinedList[iter].col,
                    combinedList[iter].rotation);
        }

        //Stop timer
        after = wallClockTime();

        printf("Parallel SETL took %1.2f seconds\n",
            ((float)(after - before))/1000000000);
    }

    MPI_Finalize();

    //Clean up
    deleteList( list );

    freeSquareMatrix( curW );
    freeSquareMatrix( nextW );

    freeSquareMatrix( patterns[0] );
    freeSquareMatrix( patterns[1] );
    freeSquareMatrix( patterns[2] );
    freeSquareMatrix( patterns[3] );

    return 0;
}

/***********************************************************
  Helper functions
***********************************************************/

void die(int lineNo)
{
    fprintf(stderr, "Error at line %d. Exiting\n", lineNo);
    exit(1);
}

long long wallClockTime( )
{
#ifdef __linux__
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}

/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue )
{
    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.
    contiguous = (char*) malloc(sizeof(char) * size * size);
    if (contiguous == NULL)
        die(__LINE__);


    memset(contiguous, defaultValue, size * size );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL)
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}

void printSquareMatrix( char** matrix, int size )
{
    int i,j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix( char** matrix )
{
    if (matrix == NULL) return;

    free( matrix[0] );
}

/***********************************************************
   World  related functions
***********************************************************/

char** readWorldFromFile( char* fname, int* sizePtr )
{
    FILE* inf;

    char temp, **world;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);

    //Using the "halo" approach
    // allocated additional top + bottom rows
    // and leftmost and rightmost rows to form a boundary
    // to simplify computation of cell along edges
    world = allocateSquareMatrix( size + 2, DEAD );

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            fscanf(inf, "%c", &world[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return world;

}

int countNeighbours(char** world, int row, int col)
//Assume 1 <= row, col <= size, no check
{
    int i, j, count;

    count = 0;
    for(i = row-1; i <= row+1; i++){
        for(j = col-1; j <= col+1; j++){
            count += (world[i][j] == ALIVE );
        }
    }

    //discount the center
    count -= (world[row][col] == ALIVE);

    return count;

}

void evolveWorld(char** curWorld, char** nextWorld,
    int size, int iter, MPI_Comm cartesian_comm)
{
    int i, j, liveNeighbours;

    char* nextTopRow = (char*)malloc(sizeof(char) * currColSize);
    char* nextBotRow = (char*)malloc(sizeof(char) * currColSize);
    char* nextLeftCol = (char*)malloc(sizeof(char) * currRowSize);
    char* nextRightCol = (char*)malloc(sizeof(char) * currRowSize);
    char nextLeftTop, nextRightTop, nextLeftBot, nextRightBot;
    MPI_Request* cartRequest = (MPI_Request*)malloc(sizeof(MPI_Request) * 16);
    MPI_Status* status = (MPI_Status*)malloc(sizeof(MPI_Status) * 16);

    if (iter > 0) {
        char* topRow = (char*)malloc(sizeof(char) * currColSize);
        char* botRow = (char*)malloc(sizeof(char) * currColSize);
        char* leftCol = (char*)malloc(sizeof(char) * currRowSize);
        char* rightCol = (char*)malloc(sizeof(char) * currRowSize);
        for (i = rowStart; i <= rowEnd; i++){
            leftCol[i-rowStart] = curWorld[i][colStart];
            rightCol[i-rowStart] = curWorld[i][colEnd];
        }

        for (j = colStart; j <= colEnd; j++) {
            topRow[j-colStart] = curWorld[rowStart][j];
            botRow[j-colStart] = curWorld[rowEnd][j];
        }

        for (i = UP; i <= DOWNRIGHT; i++) {
            if (neighbours[i] < 0) {
                cartRequest[i] = MPI_REQUEST_NULL;
                cartRequest[i+8] = MPI_REQUEST_NULL;
                continue;
            }
            if (i <= RIGHT) {
                char *sendBuf, *recBuf;
                int transferSize;
                switch (i) {
                    case UP:
                        sendBuf = topRow;
                        recBuf = nextTopRow;
                        transferSize = currColSize;
                        break;
                    case DOWN:
                        sendBuf = botRow;
                        recBuf = nextBotRow;
                        transferSize = currColSize;
                        break;
                    case LEFT:
                        sendBuf = leftCol;
                        recBuf = nextLeftCol;
                        transferSize = currRowSize;
                        break;
                    case RIGHT:
                        sendBuf = rightCol;
                        recBuf = nextRightCol;
                        transferSize = currRowSize;
                        break;
                }


                MPI_Isend(sendBuf, transferSize, MPI_CHAR, neighbours[i], getTag(iter, rank), cartesian_comm, &cartRequest[i]);

                MPI_Irecv(recBuf, transferSize, MPI_CHAR, neighbours[i], getTag(iter, neighbours[i]), cartesian_comm, &cartRequest[i + 8]);
            }
            else {
                char *sendBuf, *recBuf;
                switch (i) {
                    case UPLEFT:
                        sendBuf = &curWorld[rowStart][colStart];
                        recBuf = &nextLeftTop;
                        break;
                    case UPRIGHT:
                        sendBuf = &curWorld[rowStart][colEnd];
                        recBuf = &nextRightTop;
                        break;
                    case DOWNLEFT:
                        sendBuf = &curWorld[rowEnd][colStart];
                        recBuf = &nextLeftBot;
                        break;
                    case DOWNRIGHT:
                        sendBuf = &curWorld[rowEnd][colEnd];
                        recBuf = &nextRightBot;
                        break;
                }
                MPI_Isend(sendBuf, 1, MPI_CHAR, neighbours[i], getTag(iter, rank), cartesian_comm, &cartRequest[i]);

                MPI_Irecv(recBuf, 1, MPI_CHAR, neighbours[i], getTag(iter, neighbours[i]), cartesian_comm, &cartRequest[i + 8]);
            }
        }

        free(topRow);
        free(botRow);
        free(leftCol);
        free(rightCol);
    }

    for (i = rowStart; i <= rowEnd; i++){
        for (j = colStart; j <= colEnd; j++){
            liveNeighbours = countNeighbours(curWorld, i, j);
            nextWorld[i][j] = DEAD;

            if (curWorld[i][j] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;
        }
    }

    if (iter > 0) {
        MPI_Waitall(16, cartRequest, status);

        if (neighbours[UPLEFT] >= 0) {
            curWorld[rowStart-1][colStart-1] = nextLeftTop;
        }

        if (neighbours[UPRIGHT] >= 0){
            curWorld[rowStart-1][colEnd+1] = nextRightTop;
        }

        if (neighbours[DOWNLEFT] >= 0){
            curWorld[rowEnd+1][colStart-1] = nextLeftBot;
        }

        if (neighbours[DOWNRIGHT] >= 0){
            curWorld[rowEnd+1][colEnd+1] = nextRightBot;
        }

        if (neighbours[UP] >= 0) {
            for(i = colStart; i <= colEnd; i++) {
                curWorld[rowStart-1][i] = nextTopRow[i - colStart];
            }
        }

        if (neighbours[DOWN] >= 0) {
            for(i = colStart; i <= colEnd; i++) {
                curWorld[rowEnd+1][i] = nextBotRow[i - colStart];
            }
        }

        if (neighbours[LEFT] >= 0) {
            for(i = rowStart; i <= rowEnd; i++) {
                curWorld[i][colStart-1] = nextLeftCol[i - rowStart];
            }
        }

        if (neighbours[RIGHT] >= 0) {
            for(i = rowStart; i <= rowEnd; i++) {
                curWorld[i][colEnd+1] = nextRightCol[i - rowStart];
            }
        }

        for(i = rowStart; i <= rowEnd; i++) {
            liveNeighbours = countNeighbours(curWorld, i, colStart);
            nextWorld[i][colStart] = DEAD;

            if (curWorld[i][colStart] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][colStart] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][colStart] = ALIVE;

            liveNeighbours = countNeighbours(curWorld, i, colEnd);
            nextWorld[i][colEnd] = DEAD;

            if (curWorld[i][colEnd] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][colEnd] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][colEnd] = ALIVE;
        }

        for(i = colStart; i <= colEnd; i++) {
            liveNeighbours = countNeighbours(curWorld, rowStart, i);
            nextWorld[rowStart][i] = DEAD;

            if (curWorld[rowStart][i] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[rowStart][i] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[rowStart][i] = ALIVE;

            liveNeighbours = countNeighbours(curWorld, rowEnd, i);
            nextWorld[rowEnd][i] = DEAD;

            if (curWorld[rowEnd][i] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[rowEnd][i] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[rowEnd][i] = ALIVE;
        }
    }

    free(nextTopRow);
    free(nextBotRow);
    free(nextLeftCol);
    free(nextRightCol);
    free(cartRequest);
    free(status);
}

/***********************************************************
   Parallel related functions
***********************************************************/

int getCorrespondingEvolverOrFinderRank(int rank) {
    if (rank < numtasks/2) {
        return rank + evolverProcessors;
    }
    else {
        return rank - finderProcessors;
    }
}

/*
 * Tag is unique. It is generated from its iteration and rank.
 */
int getTag(int iter, int senderRank) {
    return iter * rowCount * colCount + rowCount * (senderRank / colCount) + (senderRank % colCount);
}

int getDiagonalCartesianNeighbour(int dir) {
    switch (dir) {
        case UPLEFT:
            if (coordinate[0] > 0 && coordinate[1] > 0) {
                return cartRank - colCount - 1;
            }
            else {
                return -1;
            }
            break;
        case UPRIGHT:
            if (coordinate[0] > 0 && coordinate[1] < colCount - 1) {
                return cartRank - colCount + 1;
            }
            else {
                return -1;
            }
            break;
        case DOWNLEFT:
            if (coordinate[0] < rowCount - 1 && coordinate[1] > 0) {
                return cartRank + colCount - 1;
            }
            else {
                return -1;
            }
            break;
        case DOWNRIGHT:
            if (coordinate[0] < rowCount - 1 && coordinate[1] < colCount - 1) {
                return cartRank + colCount + 1;
            }
            else {
                return -1;
            }
            break;
    }
    return -1;
}

void sendWorldToFinder(char ** curW, int size, int iter) {
    int i, j;
    MPI_Request req[4];
    MPI_Isend(&(curW[0][0]), (size+2)*(size+2), MPI_CHAR, getCorrespondingEvolverOrFinderRank(rank), iter, MPI_COMM_WORLD, &req[0]);

    if (neighbours[LEFT] >= 0) {
        MPI_Isend(&(curW[0][0]), (size+2)*(size+2), MPI_CHAR, getCorrespondingEvolverOrFinderRank(neighbours[LEFT]), iter, MPI_COMM_WORLD, &req[1]);
    }

    if (neighbours[UP] >= 0) {
        MPI_Isend(&(curW[0][0]), (size+2)*(size+2), MPI_CHAR, getCorrespondingEvolverOrFinderRank(neighbours[UP]), iter, MPI_COMM_WORLD, &req[2]);
    }

    if (neighbours[UPLEFT] >= 0) {
        MPI_Isend(&(curW[0][0]), (size+2)*(size+2), MPI_CHAR, getCorrespondingEvolverOrFinderRank(neighbours[UPLEFT]), iter, MPI_COMM_WORLD, &req[3]);
    }
}

void getWorldFromEvolver(char** curW, int size, int iter, int patternSize) {
    int i, j;
    MPI_Request req[4];
    MPI_Irecv(&(curW[0][0]), (size+2)*(size+2), MPI_CHAR, getCorrespondingEvolverOrFinderRank(rank), iter, MPI_COMM_WORLD, &req[0]);

    char **rightWorld = allocateSquareMatrix( size + 2, DEAD ), **botWorld = allocateSquareMatrix( size + 2, DEAD ), **rightBotWorld = allocateSquareMatrix( size + 2, DEAD );

    if (neighbours[RIGHT] >= 0) {
        MPI_Irecv(&(rightWorld[0][0]), (size+2)*(size+2), MPI_CHAR, neighbours[RIGHT], iter, MPI_COMM_WORLD, &req[1]);
    }
    else {
        req[1] = MPI_REQUEST_NULL;
    }

    if (neighbours[DOWN] >= 0) {
        MPI_Irecv(&(botWorld[0][0]), (size+2)*(size+2), MPI_CHAR, neighbours[DOWN], iter, MPI_COMM_WORLD, &req[2]);
    }
    else {
        req[2] = MPI_REQUEST_NULL;
    }

    if (neighbours[DOWNRIGHT] >= 0) {
        MPI_Irecv(&(rightBotWorld[0][0]), (size+2)*(size+2), MPI_CHAR, neighbours[DOWNRIGHT], iter, MPI_COMM_WORLD, &req[3]);
    }
    else {
        req[3] = MPI_REQUEST_NULL;
    }

    MPI_Status stats[4];
    MPI_Waitall(4, req, stats);

    if (neighbours[RIGHT] >= 0) {
        int nextRowStart = rowStart;
        int nextRowEnd = rowEnd;
        int nextColStart = colEnd + 1;
        for (j = nextRowStart; j <= nextRowEnd; j++) {
            for (i = nextColStart; i <= nextColStart+patternSize; i++) {
                curW[j][i] = rightWorld[j][i];
            }
        }
    }

    if (neighbours[DOWN] >= 0) {
        int nextRowStart = rowEnd + 1;
        int nextColStart = colStart;
        int nextColEnd = colEnd;
        for (j = nextRowStart; j <= nextRowStart+patternSize; j++) {
            for (i = nextColStart; i <= nextColEnd; i++) {
                curW[j][i] = botWorld[j][i];
            }
        }
    }

    if (neighbours[DOWNRIGHT] >= 0) {
        int nextRowStart = rowEnd + 1;
        int nextColStart = colEnd + 1;
        for (j = nextRowStart; j <= nextRowStart+patternSize; j++) {
            for (i = nextColStart; i <= nextColStart+patternSize; i++) {
                curW[j][i] = rightBotWorld[j][i];
            }
        }
    }


    freeSquareMatrix(rightWorld);
    freeSquareMatrix(botWorld);
    freeSquareMatrix(rightBotWorld);
}

/***********************************************************
   Search related functions
***********************************************************/

char** readPatternFromFile( char* fname, int* sizePtr )
{
    FILE* inf;

    char temp, **pattern;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);

    pattern = allocateSquareMatrix( size, DEAD );

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            fscanf(inf, "%c", &pattern[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return pattern;
}

void rotate90(char** current, char** rotated, int size)
{
    int i, j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            rotated[j][size-i-1] = current[i][j];
        }
    }
}

void searchPatterns(char** world, int wSize, int iteration,
        char** patterns[4], int pSize, MATCHLIST* list)
{
    int dir;

    for (dir = N; dir <= W; dir++){
        searchSinglePattern(world, wSize, iteration,
                patterns[dir], pSize, dir, list);
    }
}

void searchSinglePattern(char** world, int wSize, int iteration,
        char** pattern, int pSize, int rotation, MATCHLIST* list)
{
    int wRow, wCol, pRow, pCol, match;
    int rank, numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int chunkSize = (wSize-pSize+1) / numtasks;

    for (wRow = rowStart; wRow <= (coordinate[0] == rowCount - 1 ? (rowEnd-pSize+1) : rowEnd); wRow++){
        for (wCol = colStart; wCol <= (coordinate[1] == colCount - 1 ? (colEnd-pSize+1) : colEnd); wCol++){
            match = 1;
#ifdef DEBUGMORE
            printf("S:(%d, %d)\n", wRow-1, wCol-1);
#endif
            for (pRow = 0; match && pRow < pSize; pRow++){
                for (pCol = 0; match && pCol < pSize; pCol++){
                    if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
#ifdef DEBUGMORE
                        printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
                            world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
#endif
                        match = 0;
                    }
                }
            }
            if (match){
                insertEnd(list, iteration, wRow-1, wCol-1, rotation);
#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
#endif
            }
        }
    }
}

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

int compareBaseMatch(const void* a, const void* b) {
    BASEMATCH* elem1 = (BASEMATCH*)a;
    BASEMATCH* elem2 = (BASEMATCH*)b;
    if (elem1->iteration < elem2->iteration) {
        return -1;
    }
    else if (elem1->iteration > elem2->iteration) {
        return 1;
    }
    else if (elem1->rotation < elem2->rotation) {
        return -1;
    }
    else if (elem1->rotation > elem2->rotation){
        return 1;
    }
    else if (elem1->row < elem2->row) {
        return -1;
    }
    else if (elem1->row > elem2->row) {
        return 1;
    }
    else if (elem1->col < elem2->col) {
        return -1;
    }
    else {
        return 1;
    }
}

MATCHLIST* newList()
{
    MATCHLIST* list;

    list = (MATCHLIST*) malloc(sizeof(MATCHLIST));
    if (list == NULL)
        die(__LINE__);

    list->nItem = 0;
    list->tail = NULL;

    return list;
}

void deleteList( MATCHLIST* list)
{
    MATCH *cur, *next;
    int i;
    //delete items first

    if (list->nItem != 0 ){
        cur = list->tail->next;
        next = cur->next;
        for( i = 0; i < list->nItem; i++, cur = next, next = next->next ) {
            free(cur);
        }

    }
    free( list );
}

void insertEnd(MATCHLIST* list,
        int iteration, int row, int col, int rotation)
{
    MATCH* newItem;

    newItem = (MATCH*) malloc(sizeof(MATCH));
    if (newItem == NULL)
        die(__LINE__);

    newItem->iteration = iteration;
    newItem->row = row;
    newItem->col = col;
    newItem->rotation = rotation;

    if (list->nItem == 0){
        newItem->next = newItem;
        list->tail = newItem;
    } else {
        newItem->next = list->tail->next;
        list->tail->next = newItem;
        list->tail = newItem;
    }

    (list->nItem)++;

}

void printList(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    printf("List size = %d\n", list->nItem);


    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        printf("%d:%d:%d:%d\n",
                cur->iteration, cur->row, cur->col, cur->rotation);
    }
}
