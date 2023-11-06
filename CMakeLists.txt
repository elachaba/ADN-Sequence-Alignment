cmake_minimum_required(VERSION 3.0)

project(distanceEdition)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS "-g -Wall")

set(REFDIR ${CMAKE_SOURCE_DIR})
set(SRCDIR ${REFDIR}/src)
set(BINDIR ${REFDIR}/bin)
set(DOCDIR ${REFDIR}/doc)
set(TESTDIR ${REFDIR}/tests)
set(REPORTDIR ${REFDIR}/report)

set(CURRTEST "Needleman-Wunsch-CA.o")

add_executable(distanceEdition ${BINDIR}/distanceEdition ${BINDIR}/${CURRTEST} ${SRCDIR}/distanceEdition.c)

add_executable(distanceEdition-perf ${BINDIR}/distanceEdition-perf ${BINDIR}/${CURRTEST} ${SRCDIR}/distanceEdition.c)
target_compile_definitions(distanceEdition-perf PRIVATE __PERF_MESURE__)

add_executable(extract-fasta-sequences-size ${BINDIR}/extract-fasta-sequences-size ${SRCDIR}/extract-fasta-sequences-size.c)

add_library(Needleman-Wunsch-recmemo ${SRCDIR}/Needleman-Wunsch-recmemo.c)
target_include_directories(Needleman-Wunsch-recmemo PUBLIC ${SRCDIR})

add_library(Needleman-Wunsch-itermemo ${SRCDIR}/Needleman-Wunsch-itermemo.c)
target_include_directories(Needleman-Wunsch-itermemo PUBLIC ${SRCDIR})

add_library(Needleman-Wunsch-CA ${SRCDIR}/Needleman-Wunsch-CA.c)
target_include_directories(Needleman-Wunsch-CA PUBLIC ${SRCDIR})

add_library(Needleman-Wunsch-CO ${SRCDIR}/Needleman-Wunsch-CO.c)
target_include_directories(Needleman-Wunsch-CO PUBLIC ${SRCDIR})

add_executable(needleman-wunsch ${SRCDIR}/needleman-wunsch-main.c)
target_link_libraries(needleman-wunsch PRIVATE Needleman-Wunsch-recmemo Needleman-Wunsch-itermemo Needleman-Wunsch-CA Needleman-Wunsch-CO)
target_include_directories(needleman-wunsch PRIVATE ${SRCDIR})

#add_subdirectory(doc)
#add_subdirectory(tests)