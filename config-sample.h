// Copy this file to "config.h" and edit as needed

//#define USE_GMP
//#define MEMORY_POOL_INITIAL_SIZE   478416994uLL  // in bytes; this particular size is enough to reach n=15
#define MEMORY_POOL_INITIAL_SIZE    2688246558uLL  // in bytes; this particular size is enough to reach n=16
//#define MEMORY_POOL_INITIAL_SIZE 14378982942uLL  // in bytes; this particular size is enough to reach n=17
#define MEMORY_POOL_GROW_RATIO 1/16  // what proportion of the memory size to grow it by when more space is needed
#define HASH_TABLE_RATIO 6
//#define SHOW_PROGRESS 16  // if defined, show progress starting at this term
//#define PRINT_POLYTETS
//#define PRINT_POLYTETS_EDGE_CASES    // print edge cases which require the latest fix to the overlap-checking algorithm
//#define PRINT_POLYTETS_WITH_SYMMETRY
//#define PRINT_SYMMETRY_TOTALS
//#define DISABLE_OVERLAP_CHECKING

#define MAXIMUM_TETCOUNT 17 // 28

//#define SORT_POLYTETS
//#define WRITE_TO_FILES
//#define RESUME_FROM_FILE

#define FILE_CHUNK_SIZE (1uLL << 30)  // needs to be less than 1<<31

//#define ALTERNATIVE_OVERLAP_CHECKER  // use the overlap checker based on Level River St's clashcheck()

// MULTITHREADING will enable threading and synchronization code. Note that with this enabled, MEMORY_POOL_INITIAL_SIZE needs to be fully preallocated (it won't grow on demand).
//#define MULTITHREADING
#define WORKER_THREADS 16
#define MAXIMUM_WORK_ASSIGNMENT 512
#define HASH_TABLE_SHARDS 512
//#define SYNC_SPINLOCK  // may be faster in the case of locks that typically last for a very short time, but uses much higher CPU when locks last longer
