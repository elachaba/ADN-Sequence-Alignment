/**
 * \file Needleman-Wunsch-CO.h
 * \brief recursive cache oblivious implementation of Needleman-Wunsch global alignment algorithm that computes the distance
 * between two genetic sequences
 * \date 28/10/2023
 * \author Aymane EL-ACHAB (Ensimag, Grenoble-INP - University Grenoble-Alpes) aymane.el-achab@grenoble-inp.org
 * \author Achraf Benabach (Ensimag, Grenoble-INP - University Grenoble-Alpes) achraf.benabach@grenoble-inp.org
 */

#include "costs.h" /* For the constant costs */

/********************************************************************************
 * Iterative Cache Aware implementation of NeedlemanWunsch
 * \fn long EditDistance_NW_CA(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A
 * \param lengthA :  number of elements in A
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B
 * \return :  edit distance between A and B }
 *
 *
 * If lengthA < lengthB, the sequences A and B are swapped.
 *
 */
long EditDistance_NW_CO(char *A, size_t lengthA, char *B, size_t lengthB);
