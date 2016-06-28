#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/remoteApi/extApi.o \
	${OBJECTDIR}/remoteApi/extApiCustom.o \
	${OBJECTDIR}/remoteApi/extApiPlatform.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-lpthread

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/markov

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/markov: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/markov ${OBJECTFILES} ${LDLIBSOPTIONS} `pkg-config --libs opencv`

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DMAX_EXT_API_CONNECTIONS=255 -DNON_MATLAB_PARSING -Iinclude -IremoteApi -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/remoteApi/extApi.o: remoteApi/extApi.c 
	${MKDIR} -p ${OBJECTDIR}/remoteApi
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DMAX_EXT_API_CONNECTIONS=255 -DNON_MATLAB_PARSING -Iinclude -IremoteApi -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/remoteApi/extApi.o remoteApi/extApi.c

${OBJECTDIR}/remoteApi/extApiCustom.o: remoteApi/extApiCustom.c 
	${MKDIR} -p ${OBJECTDIR}/remoteApi
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DMAX_EXT_API_CONNECTIONS=255 -DNON_MATLAB_PARSING -Iinclude -IremoteApi -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/remoteApi/extApiCustom.o remoteApi/extApiCustom.c

${OBJECTDIR}/remoteApi/extApiPlatform.o: remoteApi/extApiPlatform.c 
	${MKDIR} -p ${OBJECTDIR}/remoteApi
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DMAX_EXT_API_CONNECTIONS=255 -DNON_MATLAB_PARSING -Iinclude -IremoteApi -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/remoteApi/extApiPlatform.o remoteApi/extApiPlatform.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/markov

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
