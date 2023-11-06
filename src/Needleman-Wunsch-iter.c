#include "Needleman-Wunsch-iter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strchr */

#include "characters_to_base.h" /* mapping from char to base */



/**
 * computeCost - computes the cost of aligning two characters from a sequence
 * when i < M and j < N
 * @param ctx holds the information about the sequences
 * @param i  the index of the character in sequence X
 * @param j the index of the character in the sequence Y
 * @param tmp the value of phi(i+1,j+1)
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
        res = ctx.slidingCol[j];
    }
    else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
    {  ManageBaseError( Yj );
        /* phi(i, j + 1) */
        res = ctx.slidingCol[j + 1];
    }
    else
    {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */
        long min = /* initialization  with cas 1*/
                ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST
                                    : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST )
                )
                + /* ph(i + 1, j + 1) */ tmp;
        { long cas2 = INSERTION_COST + ctx.slidingCol[j] ;
            if (cas2 < min) min = cas2 ;
        }
        { long cas3 = INSERTION_COST + ctx.slidingCol[j + 1] ;
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
    {
        ctx.X = A ;
        ctx.M = lengthA ;
        ctx.Y = B ;
        ctx.N = lengthB ;
    }
    else
    {
        ctx.X = B ;
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
    long res, min, tmp; /* To help computing and storing the final result */
    size_t M, N; /* M (resp. N) is the size of X (resp. Y) */

    /* We initiate the X and Y sequences,
     * X is the longest sequence, Y the shortest
     */

    struct NW_MemoIter ctx = initSequences(A, lengthA, B, lengthB);
    M = ctx.M ;
    N = ctx.N ;

    /*  Allocation and initialization of the sliding column */
    ctx.slidingCol = malloc((N + 1) * sizeof(long));
    if (ctx.slidingCol == NULL) { perror("EditDistance_NW_Iter: malloc of ctx_memB" ); exit(EXIT_FAILURE); }

    ctx.slidingCol[N] = 0;
    for (int j = N - 1; j >= 0; j--)
    {
        Yj = ctx.Y[j];
        ctx.slidingCol[j] = (isBase(Yj) ? INSERTION_COST : 0) + ctx.slidingCol[j + 1];
    }

    /* We proceed by computing columns, and overwriting the previous one
     * only the last updated column is needed to create the next column, which costs only O(N) in term of memory.
     */
    for (int i = M - 1; i >= 0; i--)
    {
        /* Compute phi(i,N) after storing the existing value in tmp */
        tmp = ctx.slidingCol[N];
        Xi = ctx.X[i];
        ctx.slidingCol[N] = (isBase(Xi) ? INSERTION_COST : 0) + ctx.slidingCol[N];
        /* Compute the current column */
        for (int j = N - 1; j >= 0; j--)
        {
            min = computeCost(ctx, i, j, tmp);
            tmp = ctx.slidingCol[j]; /* Update the value of phi(i + 1, j + 1), i is fixed. */
            ctx.slidingCol[j] = min; /* Update the value in the column, new phi(i, j) */
        }
    }

    /* Load phi(0,0) */
    res = ctx.slidingCol[0];

    /* Deallocation of ctx.slidingCol */
    free(ctx.slidingCol);

    return res;
}
