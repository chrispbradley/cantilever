!> Main program
PROGRAM CantileverExample

  USE OpenCMISS
 
  IMPLICIT NONE

  !Test program parameters
  REAL(OC_RP), PARAMETER :: WIDTH=60.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=40.0_OC_RP
  REAL(OC_RP), PARAMETER :: HEIGHT=40.0_OC_RP
  
  REAL(OC_RP), PARAMETER :: DENSITY=9.0E-4_OC_RP !in g mm^-3
  REAL(OC_RP), PARAMETER :: GRAVITY(3)=[0.0_OC_RP,0.0_OC_RP,-9.8_OC_RP] !in m s^-2
  INTEGER(OC_Intg), PARAMETER :: NUMBER_OF_LOAD_INCREMENTS=2

  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DISPLACEMENT_BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: PRESSURE_BASIS_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_GEOMETRY_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: FIELD_FIBRE_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: FIELD_MATERIAL_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: FIELD_DEPENDENT_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: FIELD_SOURCE_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program variables
  INTEGER(OC_Intg) :: displacementInterpolationType
  INTEGER(OC_Intg) :: pressureInterpolationType
  INTEGER(OC_Intg) :: pressureMeshComponent
  INTEGER(OC_Intg) :: numberOfGaussXi
  INTEGER(OC_Intg) :: scalingType
  
  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex
  INTEGER(OC_Intg) :: numberOfComputationalNodes,numberOfDomains,computationalNodeNumber
  INTEGER(OC_Intg) :: nodeNumber,nodeDomain,nodeIdx,componentIdx,derivativeIdx,leftNormalXi
  INTEGER(OC_Intg), ALLOCATABLE :: leftSurfaceNodes(:)
  INTEGER(OC_Intg) :: numberOfArguments,argumentLength,argumentStatus
  CHARACTER(LEN=255) :: commandArgument

  !CMISS variables
  TYPE(OC_BasisType) :: displacementBasis,pressureBasis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: geometricField,fibreField,materialField,dependentField,sourceField,equationsSetField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: solver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic CMISS variables
  INTEGER(OC_Intg) :: err

  LOGICAL  :: directoryExists = .FALSE.

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  CALL OC_OutputSetOn("Cantilever",err)
  !Intialise OpenCMISS
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Read in arguments and overwrite default values
  !Usage: CantileverExample [Displacement Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  displacementInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  numberOfGlobalXElements=3
  numberOfGlobalYElements=2
  numberOfGlobalZElements=2
  scalingType=OC_FIELD_ARITHMETIC_MEAN_SCALING

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
  IF(displacementInterpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    IF(numberOfArguments >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,commandArgument,argumentLength,argumentStatus)
      IF(argumentStatus>0) CALL HandleError("Error for command argument 5.")
      READ(commandArgument(1:argumentLength),*) scalingType
      IF(scalingType<0.OR.scalingType>5) CALL HandleError("Invalid scaling type.")
    ENDIF
  ELSE
    scalingType=OC_FIELD_NO_SCALING
  ENDIF
  SELECT CASE(displacementInterpolationType)
  CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=2
    pressureMeshComponent=1
  CASE(OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=3
    pressureMeshComponent=2
    pressureInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE(OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION,OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
    numberOfGaussXi=4
    pressureMeshComponent=2
    !Should generally use quadratic interpolation but use linear to match CMISS example 5e
    pressureInterpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE DEFAULT
    numberOfGaussXi=0
    pressureMeshComponent=1
  END SELECT
  WRITE(*,'("Interpolation: ", i3)') displacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') scalingType

  !Get the number of computational nodes and this computational node number
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !Create a 3D rectangular cartesian coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !Create a region and assign the coordinate system to the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"Region",err)
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  CALL OC_Region_CreateFinish(region,err)

  !Define basis function for displacement
  CALL OC_Basis_Initialise(displacementBasis,err)
  CALL OC_Basis_CreateStart(DISPLACEMENT_BASIS_USER_NUMBER,context,displacementBasis,err)
  SELECT CASE(displacementInterpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(displacementBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(displacementBasis,OC_BASIS_SIMPLEX_TYPE,err)
  END SELECT
  CALL OC_Basis_NumberOfXiSet(displacementBasis,3,err)
  CALL OC_Basis_InterpolationXiSet(displacementBasis,[displacementInterpolationType,displacementInterpolationType, &
      & displacementInterpolationType],err)
  IF(numberOfGaussXi>0) THEN
    CALL OC_Basis_QuadratureNumberOfGaussXiSet(displacementBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
  ENDIF
  CALL OC_Basis_CreateFinish(displacementBasis,err)

  IF(pressureMeshComponent/=1) THEN
    !Basis for pressure
    CALL OC_Basis_Initialise(pressureBasis,err)
    CALL OC_Basis_CreateStart(PRESSURE_BASIS_USER_NUMBER,context,pressureBasis,err)
    SELECT CASE(pressureInterpolationType)
    CASE(1,2,3,4)
      CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
    CASE(7,8,9)
      CALL OC_Basis_TypeSet(pressureBasis,OC_BASIS_SIMPLEX_TYPE,err)
    END SELECT
    CALL OC_Basis_NumberOfXiSet(pressureBasis,3,err)
    CALL OC_Basis_InterpolationXiSet(pressureBasis,[pressureInterpolationType,pressureInterpolationType, &
        & pressureInterpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(pressureBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
    CALL OC_Basis_CreateFinish(pressureBasis,err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  IF(pressureMeshComponent==1) THEN
    CALL OC_GeneratedMesh_BasisSet(generatedMesh,[displacementBasis],err)
  ELSE
    CALL OC_GeneratedMesh_BasisSet(generatedMesh,[displacementBasis,pressureBasis],err)
  ENDIF
  !Define the mesh on the region
  IF(numberOfGlobalXElements==0) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,LENGTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  CALL OC_Decomposition_CreateFinish(decomposition,err)
 
  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !Create a field to put the geometry (defualt is geometry)
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(FIELD_GEOMETRY_USER_NUMBER,region,geometricField,err)
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  CALL OC_Field_VariableLabelSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,"Geometry",err)
  CALL OC_Field_ScalingTypeSet(geometricField,scalingType,err)
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !Create a fibre field and attach it to the geometric field
  CALL OC_Field_Initialise(fibreField,err)
  CALL OC_Field_CreateStart(FIELD_FIBRE_USER_NUMBER,region,fibreField,err)
  CALL OC_Field_TypeSet(fibreField,OC_FIELD_FIBRE_TYPE,err)
  CALL OC_Field_DecompositionSet(fibreField,decomposition,err)
  CALL OC_Field_GeometricFieldSet(fibreField,geometricField,err)
  CALL OC_Field_VariableLabelSet(fibreField,OC_FIELD_U_VARIABLE_TYPE,"Fibre",err)
  CALL OC_Field_ScalingTypeSet(fibreField,scalingType,err)
  CALL OC_Field_CreateFinish(fibreField,err)

  !Create the equations_set
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,fibreField,[OC_EQUATIONS_SET_ELASTICITY_CLASS, &
    & OC_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,OC_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !Create the dependent field
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,FIELD_DEPENDENT_USER_NUMBER,dependentField,err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"Dependent",err)
  DO componentIdx=1,3
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4,pressureMeshComponent,err)
  CALL OC_Field_ComponentMeshComponentSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4,pressureMeshComponent,err)
  IF(pressureMeshComponent==1) THEN
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4, &
      & OC_FIELD_ELEMENT_BASED_INTERPOLATION,err)
  ELSE
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,4,OC_FIELD_NODE_BASED_INTERPOLATION,err)
    CALL OC_Field_ComponentInterpolationSet(dependentField,OC_FIELD_T_VARIABLE_TYPE,4, &
      & OC_FIELD_NODE_BASED_INTERPOLATION,err)
  ENDIF
  CALL OC_Field_ScalingTypeSet(dependentField,scalingType,err)
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !Create the material field
  CALL OC_Field_Initialise(materialField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,FIELD_MATERIAL_USER_NUMBER,materialField,err)
  CALL OC_Field_VariableLabelSet(materialField,OC_FIELD_U_VARIABLE_TYPE,"Material",err)
  CALL OC_Field_VariableLabelSet(materialField,OC_FIELD_V_VARIABLE_TYPE,"Density",err)
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,2.0_OC_RP,err)
  CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,6.0_OC_RP,err)
  CALL OC_Field_ComponentValuesInitialise(materialField,OC_FIELD_V_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,DENSITY,err)

  !Create the source field with the gravity vector
  CALL OC_Field_Initialise(sourceField,err)
  CALL OC_EquationsSet_SourceCreateStart(equationsSet,FIELD_SOURCE_USER_NUMBER,sourceField,err)
  CALL OC_Field_ScalingTypeSet(sourceField,scalingType,err)
  CALL OC_EquationsSet_SourceCreateFinish(equationsSet,err)
  DO componentIdx=1,3
    CALL OC_Field_ComponentValuesInitialise(sourceField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
        & componentIdx,GRAVITY(componentIdx),err)
  ENDDO

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,2,err)
  CALL OC_Field_ParametersToFieldParametersComponentCopy(geometricField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,3,err)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,4, &
    & 0.0_OC_RP,err)
  CALL OC_Field_ParameterSetUpdateStart(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)
  CALL OC_Field_ParameterSetUpdateFinish(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,err)

  !Define the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_ELASTICITY_CLASS,OC_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & OC_PROBLEM_STATIC_FINITE_ELASTICITY_SUBTYPE],problem,err)
  CALL OC_Problem_CreateFinish(problem,err)

  !Create the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  CALL OC_ControlLoop_MaximumIterationsSet(controlLoop,NUMBER_OF_LOAD_INCREMENTS,err)
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !Create the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  CALL OC_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_NewtonSolutionToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_NewtonRelativeToleranceSet(solver,1.0E-14_OC_RP,err)
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !Create the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)

  CALL OC_GeneratedMesh_SurfaceGet(generatedMesh,OC_GENERATED_MESH_REGULAR_LEFT_SURFACE,leftSurfaceNodes,leftNormalXi,err)

  !Fix x=0 nodes in x, y and z
  DO nodeIdx=1,SIZE(leftSurfaceNodes,1)
    nodeNumber=leftSurfaceNodes(nodeIdx)
    CALL OC_Decomposition_NodeDomainGet(decomposition,1,nodeNumber,nodeDomain,err)
    IF(nodeDomain==computationalNodeNumber) THEN
      DO componentIdx=1,3
        CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber, &
          & componentIdx,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        IF(displacementInterpolationType==OC_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
          DO derivativeIdx=OC_GLOBAL_DERIV_S2,OC_GLOBAL_DERIV_S1_S2_S3
            CALL OC_BoundaryConditions_AddNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,derivativeIdx, &
              & nodeNumber,componentIdx,OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !Solve problem
  CALL OC_Problem_Solve(problem,err)

  INQUIRE(file="./results", EXIST=directoryExists)
  IF (.NOT.directoryExists) THEN
    CALL execute_command_line ("mkdir ./results")
  END IF

  !Output solution
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,"./results/Cantilever","FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,"./results/Cantilever","FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    
    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
    
  END SUBROUTINE HandleError

END PROGRAM CantileverExample

