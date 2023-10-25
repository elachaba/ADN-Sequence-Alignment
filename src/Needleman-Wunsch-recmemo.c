/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;

/** \struct NW_MemoIter
 * \brief data for memoization of recursive Needleman-Wunsch algorithm
*/
struct NW_MemoIter
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long *memoA;
    long *memoB;
} ;


/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }
   }    
   
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
    
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

/**
 * computeCost - computes the cost of aligning two characters from a sequence
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param i  the index of the character in sequence X
 * @param j the index of the character in the sequence Y
 * @return The cost of aligning Xi and Yj
 */

long computeCost(struct NW_MemoIter ctx, int i, int j)
{
    char Xi, Yj;
    long res;

    Xi = ctx.X[i];
    Yj = ctx.Y[j];
    if (!isBase(Xi))  /* skip character in Xi that is not a base */
    {
        ManageBaseError( Xi );
        /* phi(i + 1, j) */
        res = ctx.memoB[j];
    }
    else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
    {  ManageBaseError( Yj );
        /* phi(i, j + 1) */
        res = ctx.memoB[j + 1];
    }
    else
    {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */
        long min = /* initialization  with cas 1*/
                ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST
                                    : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST )
                )
                + /* ph(i + 1, j + 1) */ ctx.memoA[i + 1];
        { long cas2 = INSERTION_COST + ctx.memoB[j] ;
            if (cas2 < min) min = cas2 ;
        }
        { long cas3 = INSERTION_COST + ctx.memoB[j + 1] ;
            if (cas3 < min) min = cas3 ;
        }
        res = min ;
    }

    return (res);

}

/**
 * initSequence - initialises the the sequences X and Y
 * The sequence X is always longer the the sequence Y
 * @param A The sequence from the first file
 * @param lengthA
 * @param B The sequence from the second file
 * @param lengthB
 * @return NW_MemoIter where we store all the info about the
 * sequences
 */

struct NW_MemoIter initSequences(char *A, size_t lengthA, char *B, size_t lengthB)
{
    struct NW_MemoIter ctx;
    if (lengthA >= lengthB)
    {  ctx.X = A ;
        ctx.M = lengthA ;
        ctx.Y = B ;
        ctx.N = lengthB ;
    }
    else
    {  ctx.X = B ;
        ctx.M = lengthB ;
        ctx.Y = A ;
        ctx.N = lengthA ;
    }
    return (ctx);
}

/**
 * EditDistance_NW_Iter:  is the main function to call, cf .h for specification
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the
 * recursivefunction EditDistance_NW_RecMemo
 * See .h file for documentation
 * @param A The sequence from the first file
 * @param lengthA
 * @param B the sequence from the second file
 * @param lengthB
 * @return return the Edition Distance between the sequences A and B
 */

long EditDistance_NW_Iter(char *A, size_t lengthA, char *B, size_t lengthB)
{
    _init_base_match() ;

    /* We define the variables to use in the function */
    char Xi, Yj; /* To store characters from the sequence x (resp. Y) */
    long res; /* To store final result */
    size_t M, N; /* M (resp. N) is the size of X (resp. Y) */

    /* We initiate the X and Y sequences,
     * X is the longest sequence, Y the shortest
     */

    struct NW_MemoIter ctx = initSequences(A, lengthA, B, lengthB);

    M = ctx.M ;
    N = ctx.N ;

    /* memoA will be used to store the lines calculated from the phi matrix */
    ctx.memoA = malloc((M + 1) * sizeof(long));
    if (ctx.memoA == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memoA" ); exit(EXIT_FAILURE); }
    /* memoB will be used to store the first column by setting i to M */
    ctx.memoB = malloc((N + 1) * sizeof(long));
    if (ctx.memoB == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memB" ); exit(EXIT_FAILURE);}

    ctx.memoA[M] = 0;
    ctx.memoB[N] = 0;
    for (int i = M - 1; i >= 0; i--)
    {
        // Pour j = N
        Xi = ctx.X[i];
        ctx.memoA[i] = (isBase(Xi) ? INSERTION_COST: 0) + ctx.memoA[i + 1];
    }

    for (int j = N - 1; j >= 0; j--)
    {
        Yj = B[j];
        ctx.memoB[j] = (isBase(Yj) ? INSERTION_COST: 0) + ctx.memoB[j + 1];
    }

    /* We proceed by computing columns, and storing the new column in the last one
     * only the last updated column is needed to create the next column, which costs only O(m).
     */
    for (int i = M - 1; i >= 0; i--)
    {
        ctx.memoB[N] = ctx.memoA[i];
        for (int j = N - 1; j >= 0; j--)
        {
            long min = computeCost(ctx, i, j);
            ctx.memoA[i + 1] = ctx.memoB[j]; /* Update the value of phi(i + 1, j + 1), i is fixed. */
            ctx.memoB[j] = min; /* Update the value in the column, new phi(i, j) */
        }
    }

    res = ctx.memoB[0];
    free(ctx.memoA);
    free(ctx.memoB);

    return res;
}

