#include "Needleman-Wunsch-CO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>/* for strchr */
#include <stdbool.h>
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/* The threshold to stop the recursion */
#define S 200


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
 * computeBlock: Computes the value of a block of the ph(i, j) matrix
 * @param ctx A variable storing the information about the sequences
 * @param tmp The last block line computed
 * @param begin_i the starting index in the sequence X
 * @param end_i the ending index in the sequence X
 * @param begin_j the starting index in the sequence Y
 * @param end_j the ending index in the sequence Y
 */

void computeBlock(struct NW_MemoIter ctx, long *tmp, int begin_i, int end_i, int begin_j, int end_j)
{
    long min;

    for (int i = begin_i; i > end_i; i--)
    {
        for (int j = begin_j; j > end_j; j--)
        {
            if (j == begin_j)
                /* When computing the values of the first line of the block
                 * we use the values of ph(i, j + 1) stored in tmp */
                min = computeCost(ctx, i, j, tmp[i], true);
            else
                min = computeCost(ctx, i, j, 0, false);
            /* we store ph(i + 1, j) to use it in the next iteration as ph(i + 1, j + 1)*/
            ctx.memoLine[i + 1] = ctx.slidingCol[j];
            ctx.slidingCol[j] = min;
        }
        tmp[i] = min; //Store the last value of the column in tmp line
    }
}

/**
 * EditDistance_NW_Rec_CO: calculates recursively the values of the matrix
 * of phi(i, j)
 * @param ctx A variable storing the information about the sequences
 * @param begin_i the starting index in the sequence X
 * @param end_i the ending index in the sequence X
 * @param begin_j the starting index in the sequence Y
 * @param end_j the ending sequence index in the sequence Y
 */

void EditDistance_NW_Rec_CO(struct NW_MemoIter ctx, long *tmp, int begin_i, int end_i, int begin_j, int end_j)
{
    int ni = begin_i - end_i, nj = begin_j - end_j;

    // If we reach S, we run the iterative program
    if ((ni <= S) && (nj <= S))
    {
        end_i = (end_i == 0 ? -1 : end_i);
        end_j = (end_j == 0 ? -1 : end_j);
        computeBlock(ctx, tmp, begin_i, end_i, begin_j, end_j);
    }
    else
    {
        /* We split the portion to compute by the biggest dimesion*/
        if (ni > nj)
        {
            EditDistance_NW_Rec_CO(ctx, tmp, begin_i, (begin_i + end_i) / 2, begin_j, end_j);
            EditDistance_NW_Rec_CO(ctx, tmp, (begin_i + end_i) / 2, end_i, begin_j, end_j);
        }
        else
        {
            EditDistance_NW_Rec_CO(ctx, tmp, begin_i, end_i, begin_j, (begin_j + end_j) / 2);
            EditDistance_NW_Rec_CO(ctx, tmp, begin_i, end_i, (begin_j  + end_j)/ 2, end_j);
        }
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

long EditDistance_NW_CO(char *A, size_t lengthA, char *B, size_t lengthB)
{
    _init_base_match();

    char Xi, Yj;
    long res, *tmp;
    size_t M, N;
    struct NW_MemoIter ctx = initSequences(A, lengthA, B, lengthB);

    M = ctx.M;
    N = ctx.N;

    /* Allocating the memory for memoLine and tmp and initialising them */
    ctx.memoLine = malloc((M + 1) * sizeof(long));
    if (ctx.memoLine == NULL) { perror("EditDistance_NW_CO: malloc of memoLine failed!\n"); exit(EXIT_FAILURE); }

    tmp = malloc((M + 1) * sizeof(long));
    if (tmp == NULL) { perror("EditDistance_NW_CO: malloc of tmp\n"); exit(EXIT_FAILURE); }

    ctx.memoLine[M] = 0;
    for (int i = M - 1; i >= 0; i--)
    {
        Xi = ctx.X[i];
        ctx.memoLine[i] = (isBase(Xi) ? INSERTION_COST : 0) + ctx.memoLine[i + 1];
        tmp[i] = ctx.memoLine[i];
    }

    /* Allocating the memory for slidingCol and initialising */
    ctx.slidingCol = malloc((N + 1) * sizeof(long));
    if (ctx.memoLine == NULL) { perror("EditDistance_NW_CO: malloc of slidingCol failed!\n"); exit(EXIT_FAILURE); }

    ctx.slidingCol[N] = 0;
    for (int j = N - 1; j >= 0; j--)
    {
        Yj = ctx.Y[j];
        ctx.slidingCol[j] = ((isBase(Yj) ? INSERTION_COST : 0)) + ctx.slidingCol[j + 1];
    }

    EditDistance_NW_Rec_CO(ctx, tmp, M - 1, 0, N - 1, 0);
    /* phi(0, 0)*/
    res = ctx.slidingCol[0];

    free(ctx.slidingCol);
    free(ctx.memoLine);

    return (res);
}