
#include "Needleman-Wunsch-itermemo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */



/**
 * computeCost - computes the cost of aligning two characters from a sequence
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param i  the index of the character in sequence X
 * @param j the index of the character in the sequence Y
 * @return The cost of aligning Xi and Yj
 */

long computeCost(struct NW_MemoIter ctx, int i, int j, long tmp)
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
                + /* ph(i + 1, j + 1) */ /*ctx.memoA[i + 1]*/ tmp;
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
    long res, min, tmp; /* To store final result */
    size_t M, N; /* M (resp. N) is the size of X (resp. Y) */

    /* We initiate the X and Y sequences,
     * X is the longest sequence, Y the shortest
     */

    struct NW_MemoIter ctx = initSequences(A, lengthA, B, lengthB);

    M = ctx.M ;
    N = ctx.N ;

    /* memoA will be used to store the lines calculated from the phi matrix */

    ctx.memoB = malloc((N + 1) * sizeof(long));
    if (ctx.memoB == NULL) { perror("EditDistance_NW_Iter: malloc of ctx_memB" ); exit(EXIT_FAILURE); }

    ctx.memoB[N] = 0;

    for (int j = N - 1; j >= 0; j--)
    {

        Yj = ctx.Y[j];
        ctx.memoB[j] = (isBase(Yj) ? INSERTION_COST : 0) + ctx.memoB[j + 1];
    }

    /* We proceed by computing columns, and storing the new column in the last one
     * only the last updated column is needed to create the next column, which costs only O(n).
     */
    for (int i = M - 1; i >= 0; i--)
    {
        tmp = ctx.memoB[N];
        /*ctx.memoB[N] = ctx.memoA[i];*/
        Xi = ctx.X[i];
        ctx.memoB[N] = (isBase(Xi) ? INSERTION_COST : 0) + ctx.memoB[N];
        for (int j = N - 1; j >= 0; j--)
        {
            min = computeCost(ctx, i, j, tmp);
            tmp = ctx.memoB[j]; /* Update the value of phi(i + 1, j + 1), i is fixed. */
            ctx.memoB[j] = min; /* Update the value in the column, new phi(i, j) */
        }
    }

    res = ctx.memoB[0];
    free(ctx.memoB);

    return res;
}
