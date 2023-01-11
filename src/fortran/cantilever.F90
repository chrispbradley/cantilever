!> Main program
PROGRAM CantileverExample

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif
  
  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Test program parameters
  REAL(CMISSRP), PARAMETER :: WIDTH=60.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=40.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=40.0_CMISSRP
  
  REAL(CMISSRP), PARAMETER :: DENSITY=9.0E-4_CMISSRP !in g mm^-3
  REAL(CMISSRP), PARAMETER :: GRAVITY(3)=[0.0_CMISSRP,0.0_CMISSRP,-9.8_CMISSRP] !in m s^-2
  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_LOAD_INCREMENTS=2

  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DISPLACEMENT_BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: FIELD_SOURCE_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program variables
  INTEGER(CMISSIntg) :: displacementInterpolationType
  INTEGER(CMISSIntg) :: pressureInterpolationType
  INTEGER(CMISSIntg) :: pressureMeshComponent
  INTEGER(CMISSIntg) :: numberOfGaussXi
  INTEGER(CMISSIntg) :: scalingType
  
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(CMISSIntg) :: decompositionIndex,equationsSetIndex
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,numberOfDomains,computationalNodeNumber
  INTEGER(CMISSIntg) :: nodeNumber,nodeDomain,nodeIdx,componentIdx,derivativeIdx,leftNormalXi
  INTEGER(CMISSIntg), ALLOCATABLE :: leftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,argumentStatus
  CHARACTER(LEN=255) :: commandArgument

  !CMISS variables
  TYPE(cmfe_BasisType) :: displacementBasis,pressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_ControlLoopType) :: controlLoop
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: geometricField,fibreField,materialField,dependentField,sourceField,equationsSetField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: solver,linearSolver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: err

  LOGICAL  :: directoryExists = .FALSE.

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  CALL cmfe_OutputSetOn("Cantilever",err)
  !Intialise OpenCMISS
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  !Read in arguments and overwrite default values
  !Usage: CantileverExample [Displacement Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  displacementInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  numberOfGlobalXElements=3
  numberOfGlobalYElements=2
  numberOfGlobalZElements=2
  scalingType=CMFE_FIELD_ARITHMETIC_MEAN_SCALING

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,argumentStatus)
    IF(argumentStatus>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) displacementInterpolationType
  ENDIF
  IF(numberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,argumentStatus)
    IF(argumentStatus>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<1) CALL HandleError("Invalid number of X elements.")
  ENDIF
  IF(numberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,argumentStatus)
    IF(argumentStatus>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalYElements
    IF(numberOfGlobalYElements<1) CALL HandleError("Invalid number of Y elements.")
  ENDIF
  IF(numberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,argumentStatus)
    IF(argumentStatus>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<1) CALL HandleError("Invalid number of Z elements.")
  ENDIF
  IF(displacementInterpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    IF(numberOfArguments >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,commandArgument,argumentLength,argumentStatus)
      IF(argumentStatus>0) CALL HandleError("Error for command argument 5.")
      READ(commandArgument(1:argumentLength),*) scalingType
      IF(scalingType<0.OR.scalingType>5) CALL HandleError("Invalid scaling type.")
    ENDIF
  ELSE
    scalingType=CMFE_FIELD_NO_SCALING
  ENDIF
  SELECT CASE(displacementInterpolationType)
  CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=2
    pressureMeshComponent=1
  CASE(CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=3
    pressureMeshComponent=2
    pressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE(CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
    numberOfGaussXi=4
    pressureMeshComponent=2
    !Should generally use quadratic interpolation but use linear to match CMISS example 5e
    pressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE DEFAULT
    numberOfGaussXi=0
    pressureMeshComponent=1
  END SELECT
  WRITE(*,'("Interpolation: ", i3)') displacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') scalingType

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"Region",err)
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL cmfe_Region_CreateFinish(region,err)

  !Define basis function for displacement
  CALL cmfe_Basis_Initialise(displacementBasis,err)
  CALL cmfe_Basis_CreateStart(DISPLACEMENT_BASIS_USER_NUMBER,context,displacementBasis,err)
  SELECT CASE(displacementInterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(displacementBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(displacementBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(displacementBasis,3,err)
  CALL cmfe_Basis_InterpolationXiSet(displacementBasis,[displacementInterpolationType,displacementInterpolationType, &
      & displacementInterpolationType],err)
  IF(numberOfGaussXi>0) THEN
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(displacementBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(displacementBasis,err)

  IF(pressureMeshComponent/=1) THEN
    !Basis for pressure
    CALL cmfe_Basis_Initialise(pressureBasis,err)
    CALL cmfe_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
    SELECT CASE(pressureInterpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
    END SELECT
    CALL cmfe_Basis_NumberOfXiSet(pressureBasis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(pressureBasis,[pressureInterpolationType,pressureInterpolationType, &
        & pressureInterpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
    CALL cmfe_Basis_CreateFinish(pressureBasis,err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  IF(pressureMeshComponent==1) THEN
    CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[displacementBasis],err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[displacementBasis,pressureBasis],err)
  ENDIF
  !Define the mesh on the region
  IF(numberOfGlobalXElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,LENGTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)
 
  !Decompose
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL cmfe_Field_ScalingTypeSet(geometricField,scalingType,err)
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(fibreField,err)
  CALL cmfe_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,err)
  CALL cmfe_Field_DecompositionSet(fibreField,decomposition,err)
  CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL cmfe_Field_ScalingTypeSet(fibreField,scalingType,err)
  CALL cmfe_Field_CreateFinish(fibreField,err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(equationsSetField,err)
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  DO componentIdx=1,3
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,pressureMeshComponent,err)
  CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,pressureMeshComponent,err)
  IF(pressureMeshComponent==1) THEN
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  ELSE
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,err)
  ENDIF
  CALL cmfe_Field_ScalingTypeSet(dependentField,scalingType,err)
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the material field
  CALL cmfe_Field_Initialise(materialField,err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,FIELD_MATERIAL_USER_NUMBER,materialField,err)
  CALL cmfe_Field_VariableLabelSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL cmfe_Field_VariableLabelSet(materialField,CMFE_FIELD_V_VARIABLE_TYPE,"Density",err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP,err)
  CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,DENSITY,err)

  !Create the source field with the gravity vector
  CALL cmfe_Field_Initialise(sourceField,err)
  CALL cmfe_EquationsSet_SourceCreateStart(equationsSet,FIELD_SOURCE_USER_NUMBER,sourceField,err)
  CALL cmfe_Field_ScalingTypeSet(sourceField,scalingType,err)
  CALL cmfe_EquationsSet_SourceCreateFinish(equationsSet,err)
  DO componentIdx=1,3
    CALL cmfe_Field_ComponentValuesInitialise(sourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & componentIdx,GRAVITY(componentIdx),err)
  ENDDO

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,err)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & 0.0_CMISSRP,err)
  CALL cmfe_Field_ParameterSetUpdateStart(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)
  CALL cmfe_Field_ParameterSetUpdateFinish(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,err)

  !Define the problem
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL cmfe_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  CALL cmfe_ControlLoop_Initialise(controlLoop,err)
  CALL cmfe_Problem_ControlLoopGet(problem,CMFE_CONTROL_LOOP_NODE,controlLoop,err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(controlLoop,NUMBER_OF_LOAD_INCREMENTS,err)
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(solver,1.0E-14_CMISSRP,err)
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  CALL cmfe_GeneratedMesh_SurfaceGet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,leftSurfaceNodes,leftNormalXi,err)

  !Fix x=0 nodes in x, y and z
  DO nodeIdx=1,SIZE(leftSurfaceNodes,1)
    nodeNumber=leftSurfaceNodes(nodeIdx)
    CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      DO componentIdx=1,3
        CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber, &
          & componentIdx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        IF(displacementInterpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
          DO derivativeIdx=CMFE_GLOBAL_DERIV_S2,CMFE_GLOBAL_DERIV_S1_S2_S3
            CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,derivativeIdx, &
              & nodeNumber,componentIdx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL cmfe_Problem_Solve(problem,err)

  INQUIRE(file="./results", EXIST=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL execute_command_line ("mkdir ./results")
  END IF

  !Output solution
  CALL cmfe_Fields_Initialise(fields,err)
  CALL cmfe_Fields_Create(region,fields,err)
  CALL cmfe_Fields_NodesExport(fields,"./results/Cantilever","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(fields,"./results/Cantilever","FORTRAN",err)
  CALL cmfe_Fields_Finalise(fields,err)

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
    
  END SUBROUTINE HandleError

END PROGRAM CantileverExample

