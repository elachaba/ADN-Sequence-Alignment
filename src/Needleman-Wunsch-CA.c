#include "Needleman-Wunsch-CA.h"
#include "characters_to_base.h" /* mapping from char to base */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


/* Z is the size of the cache */
#define Z 4096
/* K is the size of the block
 * K = sqrt(Z)
*/
#define K 64



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
 * computeCost - computes the cost of aligning two characters from a sequence
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param i  the index of the character in sequence X
 * @param j the index of the character in the sequence Y
 * @param tmp phi(i, j + 1) when useTmp is true
 * @param useTmp is true when calculating the first value of a column
 * @return The cost of aligning Xi and Yj
 */

long computeCost(struct NW_MemoIter ctx, int i, int j, long tmp, bool useTmp)
{
    char Xi, Yj;
    long res;

    Xi = ctx.X[i];
    Yj = ctx.Y[j];
    if (!isBase(Xi))  /* skip character in Xi that is not a base */
    {
        ManageBaseError( Xi );
        /* phi(i + 1, j) */
        res = ctx.slidingCol[j];
    }
    else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
    {  ManageBaseError( Yj );
        /* phi(i, j + 1) */
        res = (useTmp ? tmp : ctx.slidingCol[j + 1]);
    }
    else
    {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */
        long min = /* initialization  with cas 1*/
                ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST
                                    : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST )
                )
                + /* ph(i + 1, j + 1) */ ctx.memoLine[i + 1];
        { long cas2 = INSERTION_COST + ctx.slidingCol[j] ;
            if (cas2 < min) min = cas2 ;
        }
        { long cas3 = INSERTION_COST + (useTmp ? tmp : ctx.slidingCol[j + 1]) ;
            if (cas3 < min) min = cas3 ;
        }
        res = min ;
    }
    return (res);
}

/**
 * computeBlock - computes a block column by column
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param tmp an array to store the last line of a computed block
 * @param begin_i  the block's first column's index
 * @param end_i the index of the last column of the block
 * @param begin_j  the block's first line's index
 * @param ennd_j the index of the last line of the block
 * @return void
 */

void computeBlock(struct NW_MemoIter ctx, long *tmp, int begin_i, int end_i, int begin_j, int end_j)
{
    long min; /* The cost of aligning two elements */

    /* When working on the first line, we initiate tmp
     * with the values stored in the line N (the values
     * between column begin_i et column end_i + 1)
    */
    if (begin_j == ctx.N - 1)
    {
        for (int i = begin_i; i > end_i; i--)
        {
            tmp[K - (begin_i - i) - 1] = ctx.memoLine[i];
        }
    }
    for (int i = begin_i; i > end_i; i--)
    {
        for (int j = begin_j; j > end_j; j--)
        {
            if (j == begin_j)
                /* When computing the values of the first line of the block
                 * we use the values of ph(i, j + 1) stored in tmp */
                min = computeCost(ctx, i, j, tmp[K - (begin_i - i) - 1], true);
            else
                min = computeCost(ctx, i, j, 0, false);
            /* we store ph(i + 1, j) to use it in the next iteration as ph(i + 1, j + 1)*/
            ctx.memoLine[i + 1] = ctx.slidingCol[j];
            ctx.slidingCol[j] = min; //update the value of the column
        }
        tmp[K - (begin_i - i) - 1] = min; //Store the last value of the column in tmp line
    }
}

/**
 * EditDistance_NW_CA:  is the main function to call, cf .h for specification
 * See .h file for documentation
 * @param A The sequence from the first file
 * @param lengthA
 * @param B the sequence from the second file
 * @param lengthB
 * @return return the Edition Distance between the sequences A and B
 */

long EditDistance_NW_CA(char *A, size_t lengthA, char *B, size_t lengthB)
{
    _init_base_match();

    char Xi, Yj;
    long res, *tmp;
    size_t M, N;
    int end_i, end_j;
    struct NW_MemoIter ctx = initSequences(A, lengthA, B, lengthB);

    M = ctx.M;
    N = ctx.N;


    /*  Allocation and initialization of ctx.memoLine  */
    ctx.memoLine = malloc((M + 1) * sizeof(long));
    if (ctx.memoLine == NULL) { perror("EditDistance_NW_CA: malloc of memoLine failed!\n"); exit(EXIT_FAILURE); }
    ctx.memoLine[M] = 0;
    for (int i = M - 1; i >= 0; i--)
    {
        Xi = ctx.X[i];
        ctx.memoLine[i] = (isBase(Xi) ? INSERTION_COST : 0) + ctx.memoLine[i + 1];
    }

    /*  Allocation and initialization of the sliding column  */
    ctx.slidingCol = malloc((N + 1) * sizeof(long));
    if (ctx.slidingCol == NULL) { perror("EditDistance_NW_CA: malloc of slidingCol failed!\n"); exit(EXIT_FAILURE); }
    ctx.slidingCol[N] = 0;
    for (int j = N - 1; j >= 0; j--)
    {
        Yj = ctx.Y[j];
        ctx.slidingCol[j] = ((isBase(Yj) ? INSERTION_COST : 0)) + ctx.slidingCol[j + 1];
    }

    /*  Allocation of the tmp array  */
    tmp = malloc(K * sizeof(long));
    if (tmp == NULL) { perror("EditDistance_NW_CA: malloc of tmp\n"); exit(EXIT_FAILURE); }


    /* We proceed by computing K*K (or less) sized block, and then move to the block below */
    for (int I = M - 1; I >= 0; I -= K)
    {
        end_i = ((I - K) > 0 ? I - K : -1);
        for (int J = N - 1; J >= 0; J -= K)
        {
            end_j = ((J - K) > 0 ? J - K: -1);
            computeBlock(ctx, tmp, I, end_i, J, end_j);
        }
    }

    /* Load phi(0,0) */
    res = ctx.slidingCol[0];

    /* Deallocation */
    free(ctx.slidingCol);
    free(ctx.memoLine);
    free(tmp);

    return (res);
}