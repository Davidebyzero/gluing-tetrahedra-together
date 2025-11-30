// Copy this file to "config.h" and edit as needed

#define USE_GMP
#define MEMORY_POOL_INITIAL_SIZE (64uLL * 1024)  // in bytes; if the goal is to use more than half of available RAM, this must be preallocated at full expected size
#define MEMORY_POOL_GROW_RATIO 1/16  // what proportion of the memory size to grow it by when more space is needed
#define HASH_TABLE_RATIO 6
//#define SHOW_PROGRESS 16  // if defined, show progress starting at this term
//#define PRINT_POLYTETS               // requires USE_GMP
//#define PRINT_POLYTETS_EDGE_CASES    // requires USE_GMP; print edge cases which require the latest fix to the overlap-checking algorithm
//#define PRINT_POLYTETS_WITH_SYMMETRY // requires USE_GMP
//#define PRINT_SYMMETRY_TOTALS
//#define DISABLE_OVERLAP_CHECKING

#define MAXIMUM_TETCOUNT 17 // 28

//#define WRITE_TO_FILES
//#define RESUME_FROM_FILE

#define FILE_CHUNK_SIZE (1uLL << 30)  // needs to be less than 1<<31

// MULTITHREADING will enable threading and synchronization code. Note that with this enabled, MEMORY_POOL_INITIAL_SIZE needs to be fully preallocated (it won't grow on demand).
//#define MULTITHREADING
#define WORKER_THREADS 16
#define MAXIMUM_WORK_ASSIGNMENT 512
#define HASH_TABLE_SHARDS 512
