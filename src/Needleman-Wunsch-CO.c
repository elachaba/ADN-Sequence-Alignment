#include "Needleman-Wunsch-itermemo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>/* for strchr */
#include <stdbool.h>
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

#define S 200



/**
 * computeCost - computes the cost of aligning two characters from a sequence
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param i  the index of the character in sequence X
 * @param j the index of the character in the sequence Y
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
        res = ctx.memoB[j];
    }
    else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
    {  ManageBaseError( Yj );
        /* phi(i, j + 1) */
        res = (useTmp ? tmp : ctx.memoB[j + 1]);
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
        { long cas3 = INSERTION_COST + (useTmp ? tmp : ctx.memoB[j + 1]) ;
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
                min = computeCost(ctx, i, j, tmp[i], true);
            else
                min = computeCost(ctx, i, j, 0, false);
            ctx.memoA[i + 1] = ctx.memoB[j];
            ctx.memoB[j] = min;
        }
        tmp[i] = min;
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

    if ((ni <= S) && (nj <= S))
    {
        end_i = (end_i == 0 ? -1 : end_i);
        end_j = (end_j == 0 ? -1 : end_j);
        computeBlock(ctx, tmp, begin_i, end_i, begin_j, end_j);
    }
    else
    {
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

    ctx.memoA = malloc((M + 1) * sizeof(long));
    if (ctx.memoA == NULL) { perror("EditDistance_NW_CO: malloc of memoA failed!\n"); exit(EXIT_FAILURE); }
    ctx.memoB = malloc((N + 1) * sizeof(long));
    if (ctx.memoA == NULL) { perror("EditDistance_NW_CO: malloc of memoB failed!\n"); exit(EXIT_FAILURE); }
    tmp = malloc((M + 1) * sizeof(long));
    if (tmp == NULL) { perror("EditDistance_NW_CO: malloc of tmp\n"); exit(EXIT_FAILURE); }


    ctx.memoA[M] = 0;
    ctx.memoB[N] = 0;
    for (int i = M - 1; i >= 0; i--)
    {
        Xi = ctx.X[i];
        ctx.memoA[i] = (isBase(Xi) ? INSERTION_COST : 0) + ctx.memoA[i + 1];
        tmp[i] = ctx.memoA[i];
    }
    for (int j = N - 1; j >= 0; j--)
    {
        Yj = ctx.Y[j];
        ctx.memoB[j] = ((isBase(Yj) ? INSERTION_COST : 0)) + ctx.memoB[j + 1];
    }

    EditDistance_NW_Rec_CO(ctx, tmp, M - 1, 0, N - 1, 0);
    res = ctx.memoB[0];
    free(ctx.memoB);
    free(ctx.memoA);

    return (res);
}