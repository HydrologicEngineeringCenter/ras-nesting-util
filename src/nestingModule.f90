!=========================================================
module nestingModule
! Nesting Module
!
! Contains routines for extrapolating ADCIRC subdomains,
! subdomain water levels, and interpolating time-series,
! at points from ADCIRC domains or subdomains.
!
! Author: Alex Sanchez, IWR-HEC
!=========================================================
  implicit none
  
  double precision,parameter :: PI = acos(-1D0)
  
  type MeshType
    integer:: numElems = 0
    integer:: numNodes = 0
    integer:: numNodesParent = 0
    integer:: numElemsParent = 0
    integer,allocatable:: nodeParent(:)
    integer,allocatable:: nodeChild(:)
    integer,allocatable:: elemParent(:)
    integer,allocatable:: elemChild(:)
    integer,allocatable:: elem2node(:,:)   !Element to node connectivity
    double precision,allocatable:: nodexyz(:,:) !Global coordinates
    character(len=512):: header
    character(len=512):: filename
  endtype
  
contains

!**************************************
  subroutine deallocate_mesh(mesh)
!**************************************
    implicit none
    type(MeshType):: mesh
    
    if(allocated(mesh%nodeParent))then
      deallocate(mesh%nodeParent)
    endif
    if(allocated(mesh%nodeChild))then
      deallocate(mesh%nodeChild)
    endif
    if(allocated(mesh%elemParent))then
      deallocate(mesh%elemParent)
    endif
    if(allocated(mesh%elemChild))then
      deallocate(mesh%elemChild)
    endif
    if(allocated(mesh%elem2node))then
      deallocate(mesh%elem2node)
    endif
    if(allocated(mesh%nodexyz))then
      deallocate(mesh%nodexyz)
    endif
    
    mesh%numElems = 0
    mesh%numNodes = 0
    mesh%numNodesParent = 0
    mesh%numElemsParent = 0
    mesh%header = ''
    mesh%filename = ''
    
  endsubroutine
  
!*****************************************************************************
  subroutine nesting_gensub(parentGridFile,polyFile,childGridFile,status,msg)
!*****************************************************************************
    implicit none
    character(len=*),intent(in):: parentGridFile,polyFile,childGridFile
    integer,intent(inout):: status
    character(len=*),intent(inout):: msg
    integer:: np
    type(MeshType):: P, C
    double precision,allocatable:: polyxy(:,:)
  
    !Read Parent Mesh
    P%filename = parentGridFile
    write(*,*) 'Reading Parent Mesh: ',trim(P%filename)
    call read_grid14(P%filename, P%header, &
      P%numElems, P%numNodes, P%nodexyz, P%elem2node)
    write(*,*)'   Success'
    
    !Read Bounding Polygon
    write(*,*) 'Reading Polygon File: ',trim(polyFile)
    call read_csv(polyFile,polyxy)
    np = size(polyxy,1)
    if(np < 3)then
      msg = '  Error: Polygon must contain at least 3 points'
      status = -1
      return
    else
      write(*,*)'   Success'
    endif
    
    !Create Child Mesh
    write(*,*) 'Creating Child Mesh'
    call create_child_mesh(P, C, np, polyxy(:,1), polyxy(:,2))
    C%filename = childGridFile
    C%header   = P%header
    write(*,*)'   Success'
    
    !Write Child Mesh
    write(*,*) 'Writing Child Mesh: ',trim(C%filename)
    call write_grid14(C%filename, P%header, &
      C%numElems, C%numNodes, C%nodexyz, C%elem2node, C%nodeParent)
    write(*,*)'   Success'
    
    !Clean
    deallocate(polyxy)
    call deallocate_mesh(P)
    call deallocate_mesh(C)
    
    status = 0
    msg = 'Success'
    
  endsubroutine
  
!*****************************************************************************
  subroutine nesting_trisp(savePointFile,polyFile,childGridFile,status,msg)
!*****************************************************************************
    implicit none
    character(len=*),intent(in):: savePointFile,polyFile,childGridFile
    integer,intent(inout):: status
    character(len=*),intent(inout):: msg
    integer:: numPoints,numPoly,ierr
    double precision,allocatable:: polyxy(:,:)
    double precision,allocatable:: pointsxy(:,:)
    type(MeshType):: C
    
    !Read Save Point Coordinates
    write(*,*) 'Reading Save Point Coordinates: ',trim(savePointFile)
    call read_savepoints151(savePointFile,numPoints,pointsxy,ierr)
    write(*,*)'   Success'
    
    !Read Bounding Polygon
    write(*,*) 'Reading Polygon File: ',trim(polyFile)
    call read_csv(polyFile,polyxy)
    numPoly = size(polyxy,1)
    write(*,*)'   Success'
    
    !Triangulate Subdomain
    write(*,*) 'Creating Child Mesh'
    call triangulate_subdomain(C,numPoints,pointsxy,numPoly,polyxy,status,msg)
    write(*,*)'   Success'
    
    !Write Child Mesh
    C%filename = childGridFile
    write(*,*) 'Writing Child Mesh: ',trim(C%filename)
    call write_grid14(C%filename, 'Subdomain', &
      C%numElems, C%numNodes, C%nodexyz, C%elem2node, C%nodeParent)
    write(*,*)'   Success'
    
    !Clean
    deallocate(pointsxy,polyxy)
    call deallocate_mesh(C)
    
    status = 0
    msg = 'Success'
    
  endsubroutine
  
!***********************************************************************************
  subroutine triangulate_subdomain(C,numPoints,pointsxy,numPoly,polyxy,status,msg)
!***********************************************************************************
    use delaunay2d, only: dtris2
    implicit none
    type(MeshType),intent(inout):: C
    integer:: numPoints,numPoly
    double precision,intent(in):: pointsxy(numPoints,2)
    double precision,intent(in):: polyxy(numPoly,2)
    integer,intent(inout):: status
    character(len=*),intent(inout):: msg
    integer,allocatable:: temp(:)
    integer:: i, j, numNodes
    double precision:: xi,yi
    logical:: inPoly
    integer:: numTri
    double precision,allocatable:: nodexy(:,:)
    integer,allocatable :: tri2node(:,:),tri2tri(:,:)
    
    status = 0
    msg = 'Success'
    
    !Find Save Points within Polygon
    allocate(temp(numPoints))
    allocate(C%nodeChild(numPoints))
    C%nodeChild= 0
    numNodes = 0
    do i=1,numPoints
      xi = pointsxy(i,1)
      yi = pointsxy(i,2)
      inPoly = inPolygon(numPoly,polyxy(:,1),polyxy(:,2),xi,yi)
      if(inPoly)then
        numNodes = numNodes + 1
        temp(numNodes) = i
        C%nodeChild(i) = numNodes
      endif
    enddo
    C%numNodes = numNodes
    allocate(C%nodeParent(C%numNodes))
    C%nodeParent(1:C%numNodes) = temp(1:C%numNodes)
    deallocate(temp)
    allocate(C%nodexyz(C%numNodes,3))
    allocate(nodexy(2,numNodes))
    do i=1,C%numNodes
      j = C%nodeParent(i) !Parent Node/Save Point ID
      C%nodexyz(i,1:2) = pointsxy(j,1:2)
      C%nodexyz(i,3)   = -99999.0
      nodexy(1:2,i) = pointsxy(j,1:2)
    enddo
    
    !Triangulate
    allocate(tri2node(3,numNodes*2),tri2tri(3,numNodes*2))
    tri2node = 0
    tri2tri = 0
    call dtris2(numNodes, nodexy, numTri, tri2node, tri2tri)
    C%numElems = numTri
    allocate(C%elem2node(3,C%numElems))
    do i=1,C%numElems
      C%elem2node(1:3,i) = tri2node(1:3,i)
    enddo
    deallocate(nodexy,tri2node,tri2tri)
    
  endsubroutine
  
!************************************************
  subroutine create_child_mesh(P,C,np,xp,yp)
!************************************************
    implicit none
    type(MeshType),intent(in):: P
    type(MeshType),intent(inout):: C
    integer:: np
    double precision,intent(in):: xp(np),yp(np)
    integer,allocatable:: temp(:)
    integer:: i, j, k, numNodesChild, numElemChild
    integer:: nodes(3)
    double precision:: xi,yi
    logical:: inPoly
    
    C%numNodesParent = P%numNodes
    C%numElemsParent = P%numElems
    
    !Find Nodes within Polygon
    allocate(temp(P%numNodes))
    allocate(C%nodeChild(C%numNodesParent))
    C%nodeChild= 0
    numNodesChild = 0
    do i=1,P%numNodes
      xi = P%nodexyz(i,1)
      yi = P%nodexyz(i,2)
      inPoly = inPolygon(np,xp,yp,xi,yi)
      if(inPoly)then
        numNodesChild = numNodesChild + 1
        temp(numNodesChild) = i
        C%nodeChild(i) = numNodesChild
      endif
    enddo
    C%numNodes = numNodesChild
    allocate(C%nodeParent(C%numNodes))
    C%nodeParent(1:C%numNodes) = temp(1:C%numNodes)
    deallocate(temp)
    allocate(C%nodexyz(C%numNodes,3))
    do i=1,C%numNodes
      j = C%nodeParent(i) !Parent Node
      C%nodexyz(i,:) = P%nodexyz(j,:)
    enddo
    
    !Find Elements within Polygon
    allocate(temp(P%numElems))
    allocate(C%elemChild(P%numElems))
    C%elemChild = 0
    numElemChild = 0
    do j=1,P%numElems
      nodes = P%elem2node(:,j)
      if(all(C%nodeChild(nodes) > 0))then
        numElemChild = numElemChild + 1
        temp(numElemChild) = j
        C%elemChild(j) = numElemChild
      endif
    enddo
    C%numElems = numElemChild
    allocate(C%elemParent(C%numElems))
    C%elemParent(1:C%numElems) = temp(1:C%numElems)
    allocate(C%elem2node(3,C%numElems))
    
    do j=1,C%numElems
      k = C%elemParent(j) !Parent element
      nodes = P%elem2node(:,k) !Parent element nodes
      C%elem2node(:,j) = C%nodeChild(nodes)
    enddo
    
    !Clean
    deallocate(temp)
    
  endsubroutine
  
!**************************************************************
  subroutine nesting_extractwse(childGridFile, numWSEFiles, &
    parentWSEFiles, childWSEFiles, status, msg)
!**************************************************************
    implicit none
    character(len=*),intent(in):: childGridFile
    integer,intent(in):: numWSEFiles
    character(len=*),intent(in):: parentWSEFiles(numWSEFiles)
    character(len=*),intent(in):: childWSEFiles(numWSEFiles)
    integer,intent(inout):: status
    character(len=*),intent(inout):: msg
    type(MeshType):: C
    integer:: i
    
    !Read Mesh (Could be Parent or Child Mesh)
    C%filename = childGridFile
    write(*,*) 'Reading Mesh: ',trim(C%filename)
    call read_grid14(C%filename, C%header, C%numElems, &
      C%numNodes, C%nodexyz, C%elem2node, C%nodeParent)
    write(*,*)'   Success'
    
    !Extract Water Levels on Child Mesh
    do i=1,numWSEFiles
      write(*,*) 'Parent Water Level File: ',trim(parentWSEFiles(i))
      write(*,*) 'Child Water Level File:  ',trim(childWSEFiles(i))
      call extract_child_water_levels(parentWSEFiles(i), &
        childWSEFiles(i), C%numNodes, C%nodeParent)
    enddo
    write(*,*)'   Success'
    
    !Clean
    call deallocate_mesh(C)
    
    status = 0
    msg = 'Success'
    
  endsubroutine

!*********************************************************************************
  subroutine nesting_interpwse(childGridFile, pointsCoordFile, &
    numWSEFiles, childWSEFiles, pointsWSEFiles, factor, status, msg)
! Input:
!   childGridFile : ADCIRC Mesh (does not have to be a subdomain)
!   pointsCoordFile : Coordinates of points to be interpolated (boundary faces)
!   numWSEFiles : Number of WSE solution files
!   childWSEFiles : WSE files (do not have to be for subdomains)
!   pointsWSEFiles : Points WSE File
!*********************************************************************************
    implicit none
    character(len=*),intent(in):: childGridFile !fort.14 file
    character(len=*),intent(in):: pointsCoordFile !CSV file with format ID,X,Y or X,Y
    integer,intent(in):: numWSEFiles
    character(len=*),intent(in):: childWSEFiles(numWSEFiles) !fort.63 file
    character(len=*),intent(in):: pointsWSEFiles(numWSEFiles) !CSV file
    double precision,intent(in):: factor
    integer,intent(inout):: status
    character(len=*),intent(inout):: msg
    type(MeshType):: C
    integer:: i, npts, nCol
    double precision,allocatable:: dat(:,:)
    double precision,allocatable:: fpts(:)
    double precision,allocatable:: xpts(:),ypts(:)
    double precision:: xtrapdist = 0.0
    
    !Read Mesh (Could be Parent or Child)
    C%filename = childGridFile
    write(*,*) 'Reading Mesh: ',trim(C%filename)
    call read_grid14(C%filename, C%header, C%numElems, C%numNodes, C%nodexyz, C%elem2node)
    write(*,*)'   Success'
    
    !Read Boundary Point Coordinate File
    write(*,*) 'Reading Point Coordinate File: ',trim(pointsCoordFile)
    call read_csv(pointsCoordFile,dat)
    npts = size(dat,1)
    nCol = size(dat,2)
    allocate(fpts(npts),xpts(npts),ypts(npts))
    do i = 1,npts
      xpts(i) = dat(i,1)
      ypts(i) = dat(i,2)
    enddo
    !Compute Fraction (0-1) along boundary as headers
    fpts(1) = 0.0
    do i=2,npts
      fpts(i) = fpts(i-1) + sqrt((xpts(i) - xpts(i-1))**2 + (ypts(i) - ypts(i-1))**2)
    enddo
    fpts = fpts / fpts(npts)
    write(*,*)'   Success'
    
    write(*,*)'Number of WSE Files: ',numWSEFiles
    
    !Extract Water Levels on Points from Child Mesh
    do i=1,numWSEFiles
      write(*,*) 'Water Level Solution File: ',trim(childWSEFiles(i)) !'Child.63'
      write(*,*) 'Water Level Time Series File: ',trim(pointsWSEFiles(i)) !'Child.019'
      call interp_points_water_levels(childWSEFiles(i), pointsWSEFiles(i), &
        C%numElems, C%numNodes, C%nodexyz, C%elem2node, npts, fpts, xpts, ypts, factor, xtrapdist)
      write(*,*)'   Success'
    enddo
    
    !Clean
    if(allocated(dat))then
      deallocate(dat)
    endif
    if(allocated(fpts))then
      deallocate(fpts)
    endif
    if(allocated(xpts))then
      deallocate(xpts)
    endif
    if(allocated(ypts))then
      deallocate(ypts)
    endif
    if(allocated(dat))then
      deallocate(dat)
    endif
    call deallocate_mesh(C)
    
    status = 0
    msg = 'Success'
    
  endsubroutine

!**************************************************************************************
  subroutine interp_points_water_levels(wsefile, pointWSEFile, &
    numElems, numNodes, nodexyz, elem2node, npts, fpts, xpts, ypts, factor, xtrapdist)
! Input:
!   wsefile : Input water surface elevation file
!   pointWSEFile : Output points water surface file
!   numElems : number of elements
!   numNodes : number of nodes
!   nodexyz : node coordinates
!   elem2node : element to node connectivity
!   npts : number of points
!   ftps : Column values
!   xpts, ypts : point coordinates
!   factor : Conversion factor
!   xtrapdist : Extrapolation distance
!**************************************************************************************
    implicit none
    character(len=*),intent(in):: wsefile, pointWSEFile
    integer,intent(in):: numElems,numNodes
    double precision,intent(in):: nodexyz(numNodes,3)
    integer,intent(in):: elem2node(3,numElems)
    integer:: npts
    double precision,intent(in):: fpts(npts)
    double precision,intent(in):: xpts(npts),ypts(npts)
    double precision,intent(in):: factor, xtrapdist
    integer,allocatable:: intp(:,:)
    double precision,allocatable:: cntp(:,:)
    double precision,allocatable:: elevNodes(:), elevPts(:)
    character(len=512):: header
    integer:: i, j, node, numTimeSteps, step, ts_inc, irType
    integer:: ierr, numWetNodes, num, intTimeSteps
    double precision:: t_inc, time, elevDry, elev
    double precision,parameter:: ELEV_DRY = -99999.0
    
100 format(e22.10,i15)
150 format(f,i,i,f)
200 format(5000(e17.8,a))
250 format(5000(i17  ,a))
300 format(i11,i11,e16.7,i6,i6)
    
    !Compute Interpolation Coefficients
    allocate(intp(3,npts))
    allocate(cntp(3,npts))
    call interp_coef_tri2pts(numElems, numNodes, &
      nodexyz(:,1), nodexyz(:,2), elem2node, &
      npts, xpts, ypts, xtrapdist, intp, cntp)
  
    !Water Level File
    !RUNDES, RUNID, AGRID
    !NDSETSE, NP, DTDP*NSPOOLGE, NSPOOLGE, IRTYPE
    !TIME, IT
    !for k=1, NP
    !  k, ETA2(k)
    !end
    open(63,file=wsefile)
    read(63,'(A)') header
    read(63,*) numTimeSteps, num, t_inc, ts_inc, irType
    if(num /= numNodes)then
      write(*,*) 'ERROR: Incorrect number of nodes in water level file: ',trim(wsefile)
      stop
    endif
    
    !Point Water Level CSV File
    open(23,file=pointWSEFile,iostat=ierr)
    if(ierr /= 0)then
      write(*,*) 'ERROR: Could not open file: ',trim(pointWSEFile)
      stop
    endif
    write(23,200) ((fpts(i),','),i=1,npts-1),fpts(npts)
    
    write(*,*) 'Number of time steps: ',numTimeSteps
    write(*,*) 'Interpolating Values'
    
    intTimeSteps = numTimeSteps/20
    
    allocate(elevNodes(numNodes))
    allocate(elevPts(npts))
    do i=1,numTimeSteps
      if(numElems > 1000000 .and. mod(i,intTimeSteps)==0)then
        write(*,'(2x,I2,A1)') 100*i/numTimeSteps,'%'
      endif
      
      !Step Header
      numWetNodes = -1
      read(63,150,iostat=ierr) time, step, numWetNodes, elevDry !time in seconds
      
      !Read Node Water Levels
      if(numWetNodes > 0)then !Compact Format
        elevNodes = elevDry
        do j=1,numWetNodes
          read(63,*) node, elev
          elevNodes(node) = elev
        enddo
      else !Full format
        elevDry = ELEV_DRY
        do j=1,numNodes
          read(63,*) node,elevNodes(j)
        enddo
      endif
      
      !Interpolate Water Levels from Mesh to Points
      call interp_val_tri2pts(numNodes, elevNodes, &
        elevDry, npts, intp, cntp, elevPts)
      
      !Conversion of wet points
      if(abs(factor - 1D0) > 1d-6)then
        do j=1,npts
          if(abs(elevPts(j)-elevDry) > 1d-6)then
            elevPts(j) = elevPts(j) * factor
          endif
        enddo
      endif
      
      !Write Point Water Levels
      write(23,200) ((elevPts(j),','),j=1,npts-1),elevPts(npts)
    enddo
    
    if(numElems > 1000000)then
      write(*,'(2x,A)') '100%'
    endif
    
    close(63)
    close(23)
    
    !Clean
    deallocate(intp,cntp,elevNodes,elevPts)
    
  endsubroutine

!******************************************************
  subroutine interp_val_tri2pts(numNodes, valNodes, &
    valDry, npts, intp, cntp, valPts)
!******************************************************
    implicit none
    integer,         intent(in) :: numNodes
    double precision,intent(in) :: valNodes(numNodes)
    double precision,intent(in) :: valDry
    integer,         intent(in) :: npts
    integer,         intent(in) :: intp(3,npts)
    double precision,intent(in) :: cntp(3,npts)
    double precision,intent(out):: valpts(npts)
    integer:: i,j,idx,ind(3)
    double precision:: wcoef(3), sumw
    
    do i=1,nPts
      !Only Use Weights from Wet Nodes
      ind = intp(:,i)
      do idx=1,3
        j = ind(idx)
        if(abs(valNodes(j) - valDry) > 1e-6)then
          wcoef(idx) = cntp(idx,i) !Wet
        else
          wcoef(idx) = 0.0d0 !Dry
        endif
      enddo
      !Interpolate
      sumw = sum(wcoef)
      if(sumw > 1d-20)then !At least one node wet
        wcoef = wcoef / sumw !Normalize wet weights
        valpts(i) = sum(wcoef * valNodes(ind))
      else !All nodes wet
        valpts(i) = valDry
      endif
    enddo
    
  endsubroutine
  
!******************************************************************************************
  subroutine read_grid14(grd14file,header,numElems,numNodes,nodexyz,elem2node,nodeParent)
! Reads the parent ADCIRC grid file
!******************************************************************************************
    implicit none
    integer,    intent(out)              :: numElems       !Number of elements
    integer,    intent(out)              :: numNodes       !Number of nodes
    integer,    intent(inout),allocatable:: elem2node(:,:) !Element to node connectivity
    double precision,intent(inout),allocatable:: nodexyz(:,:)   !Node coordinates (node,coord)
    integer,         intent(inout),allocatable,optional:: nodeParent(:)
    character(len=*),intent(in)          :: grd14file      !ADCIRC grid file
    character(len=*),intent(inout)       :: header
    integer:: i,k,id,numedges,ierr,intNodes,intElems
    logical:: writeStatus = .false.
    logical:: found
  
    inquire(file=grd14file,exist=found)
    if(.not.found)then
      write(*,*) 'ERROR: Could not find ADCIRC grid file: ',trim(grd14file)
      stop
    endif
    open(unit=14,file=grd14file)
    read(14,'(A)') header
    read(14,*) numElems,numNodes
    if(numElems < numNodes)then
      write(*,'(A)') 'Error: Ivalid format. Number of elements less than number of nodes.'
      stop
    endif
    
    if(numNodes > 1000000)then
      writeStatus = .true.
      write(*,'(A)') 'Reading Nodes:'
    endif
    intNodes = numNodes/10 + 1
    intElems = numElems/10 + 1
    
    allocate(nodexyz(numNodes,3))
    if(present(nodeParent))then
      !Try first line
      allocate(nodeParent(numNodes))
      read(14,*,iostat=ierr) id,nodexyz(1,1),nodexyz(1,2),nodexyz(1,3),nodeParent(1)
      if(ierr == 0)then
        do i=2,numNodes
          if(writeStatus .and. mod(i,intNodes)==0 .and. i/=numNodes)then
            write(*,'(2x,I2,A1)') 100*i/numNodes,'%'
          endif
          read(14,*) id,nodexyz(i,1),nodexyz(i,2),nodexyz(i,3),nodeParent(i)
        enddo
      else
        nodeParent = 0 !Not a child mesh
        do i=2,numNodes
          if(writeStatus .and. mod(i,intNodes)==0 .and. i/=numNodes)then
            write(*,'(2x,I2,A1)') 100*i/numNodes,'%'
          endif
          read(14,*) id,nodexyz(i,1),nodexyz(i,2),nodexyz(i,3)
        enddo
      endif
    else
      do i=1,numNodes
        if(writeStatus .and. mod(i,intNodes)==0 .and. i/=numNodes)then
          write(*,'(2x,I2,A1)') 100*i/numNodes,'%'
        endif
        read(14,*) id,nodexyz(i,1),nodexyz(i,2),nodexyz(i,3)
      enddo
    endif
    if(writeStatus)then
      write(*,'(2x,A)') '100%'
      write(*,'(A)') 'Reading Elements:'
    endif
    
    allocate(elem2node(3,numElems))
    do i=1,numElems
      read(14,*) id,numedges,(elem2node(k,i),k=1,3)
      if(writeStatus .and. mod(i,intElems)==0 .and. i/=numElems)then
        write(*,'(2x,I2,A1)') 100*i/numElems,'%'
      endif
      if(id /= i)then
        write(*,*) 'ERROR: Problem reading ADCIRC grid connectivity at element: ',i
        stop
      endif
    enddo
    close(14)
    if(writeStatus)then
      write(*,'(2x,A)') '100%'
    endif
    
  endsubroutine
  
!*******************************************************************************************
  subroutine write_grid14(grdfile,header,numElems,numNodes,nodexyz,elem2node,nodeParent)
! Writes an ADCIRC Grid File
!*******************************************************************************************
    implicit none
    integer,         intent(in):: numElems,numNodes     !# of elements and nodes
    integer,         intent(in):: elem2node(3,numElems) !element to node connectivity
    double precision,     intent(in):: nodexyz(numNodes,3) !node coordinates (node,coord)
    integer,         intent(in),optional:: nodeParent(numNodes)
    character(len=*),intent(in):: grdfile !Grid file
    character(len=*),intent(in):: header
    integer:: i,j,k,intNodes,intElems
    logical:: writeStatus = .false.
    
222 format(I6,I6)
333 format(I6,3(1x,F20.10),1x,I10)
444 format(I6,1x,I3,3(1x,I6))
    
    if(numNodes > 1000000)then
      writeStatus = .true.
      write(*,'(A)') 'Writing Nodes:'
    endif
    intNodes = numNodes/10 + 1
    intElems = numElems/10 + 1
    
    open(144,file=grdfile)
    write(144,'(A)') trim(header)
    write(144,222) numElems,numNodes
    if(present(nodeParent))then
      do k=1,numNodes
        if(writeStatus .and. mod(i,intNodes)==0 .and. i/=numNodes)then
          write(*,'(2x,I2,A1)') 100*i/numNodes,'%'
        endif
        write(144,333) k,(nodexyz(k,i),i=1,3),nodeParent(k)
      enddo
    else
      do k=1,numNodes
        if(writeStatus .and. mod(i,intNodes)==0 .and. i/=numNodes)then
          write(*,'(2x,I2,A1)') 100*i/numNodes,'%'
        endif
        write(144,333) k,(nodexyz(k,i),i=1,3)
      enddo
    endif
    if(writeStatus)then
      write(*,'(2x,A)') '100%'
      write(*,'(A)') 'Writing Elements:'
    endif
    do j=1,numElems
      if(writeStatus .and. mod(i,intElems)==0 .and. i/=numElems)then
        write(*,'(2x,I2,A1)') 100*i/numElems,'%'
      endif
      write(144,444) j,3,(elem2node(k,j),k=1,3)
    enddo
    close(144)
    if(writeStatus)then
      write(*,'(2x,A)') '100%'
    endif
    
  endsubroutine
  
!****************************************
  subroutine read_csv(csvfile,dat)
! Reads an CSV File
!****************************************
    implicit none
    character(len=*),intent(in):: csvfile
    double precision,intent(inout),allocatable:: dat(:,:)
    integer:: i,j,ierr,numCol,numRow
    character(len=2048):: strline
    logical:: foundfile
    
    !Open File
    inquire(file=csvfile,exist=foundfile)
    if(.not.foundfile)then
      write(*,*) 'ERROR: Could not find file: ',trim(csvfile)
      pause
      stop
    endif
    open(454,file=csvfile)
  
    !Find Number of Columns
    read(454,'(A)',iostat=ierr) strline
    if(ierr/=0)then
      write(*,*) 'ERROR: Problem reading header of CSV file: ',trim(csvfile)
      pause
      stop
    endif
    numCol = count_commas(strline) + 1
    
    !Find Number of Rows
    numRow = 1
    do i=1,1000000
      read(454,*,iostat=ierr) strline
      if(ierr/=0)then
        exit
      endif
      numRow = numRow + 1
    enddo
    
    !Read Data
    rewind(454)
    allocate(dat(numRow,numCol))
    do i=1,numRow
      read(454,*,iostat=ierr) (dat(i,j),j=1,numCol)
    enddo
    
    !Close File
    close(454)
    !write(*,*) csvfile,' read succesfully'
    
  endsubroutine
  
!*******************************************************
  subroutine read_savepoints151(savePointFile,n,xy,ierr)
! Reads an ADCIRC Save Point File
!*******************************************************
    implicit none
    character(len=*),intent(in):: savePointFile
    integer,intent(out):: n
    double precision,intent(inout),allocatable:: xy(:,:)
    integer,intent(out):: ierr
    integer:: i
    logical:: foundfile
    
    ierr = 0
    
    !Open File
    inquire(file=savePointFile,exist=foundfile)
    if(.not.foundfile)then
      ierr = -1
      write(*,*) 'ERROR: Could not find file: ',trim(savePointFile)
      pause
      stop
    endif
    open(151,file=savePointFile)
    read(151,*,iostat=ierr) n
    allocate(xy(n,2))
    do i=1,n
      read(151,*,iostat=ierr) xy(i,1:2)
      if(ierr /= 0)then
        exit
      endif
    enddo
    close(151)
    
  endsubroutine
  
!************************************************
  function count_commas(string) result(num)
!************************************************
    implicit none
    integer:: i,n,num
    character(len=*),intent(inout):: string
    
    n = len_trim(string)
    num = 0
    do i=1,n
      if(string(i:i)==',')then
        num = num + 1
      endif
    enddo
  
  endfunction
  
!**********************************************************************
  subroutine interp_coef_tri2pts(numElems,numNodes,xn,yn,elem2node,&
              npts,xpts,ypts,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients from an unstructured 
! node-based triangular mesh (such as ADCIRC) to scattered points
!**********************************************************************
    implicit none
    !Parent Unstructured Triangular Mesh
    integer,         intent(in):: numElems                  !# of cells (elements)
    integer,         intent(in):: numNodes                  !# of nodes
    integer,         intent(in):: elem2node(3,numElems)     !Cells (element) to node connectivity
    double precision,intent(in):: xn(numNodes),yn(numNodes) !Node global coordinates
    !Child Points
    integer,         intent(in)   :: npts
    double precision,intent(in)   :: xpts(npts),ypts(npts)
    integer,         intent(inout):: intp(3,npts)
    double precision,intent(inout):: cntp(3,npts)
    !Extrapolation
    double precision,intent(in):: xtrapdist
    !Internal variables
    integer:: i,j,k,intPoints
    double precision:: xtri(3),ytri(3),wcoef(3),xp,yp
    double precision:: dist,distmin
    logical:: inTri
    
    !Initialize
    intp = 0
    cntp = 0d0
    
    write(*,*) 'Computing Interpolation Coefficients'
    intPoints = npts/20
    
d1: do i=1,npts !Points
      xp = xpts(i) !Global coordinates
      yp = ypts(i) !Global coordinates
      
      if(numElems > 1000000 .and. mod(i,intPoints)==0)then
        write(*,'(2x,I2,A1)') 100*i/npts,'%'
      endif
      
      !Search if point is in any triangles
      do j=1,numElems !Elements
        xtri = xn(elem2node(:,j)) !Global X-coordinates
        ytri = yn(elem2node(:,j)) !Global Y-coordinates
        inTri = inTriangle(xp,yp,xtri,ytri)
        if(inTri)then !wave point inside flow triangle
          call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
          intp(:,i) = elem2node(:,j)
          cntp(:,i) = wcoef
          cycle d1
        endif
      enddo !j
      
      !Extrapolation to nearest node
      distmin = 1d20
      do k=1,numNodes
        dist = sqrt((xp-xn(k))**2 + (yp-yn(k))**2 + 1d-15)
        if(dist < distmin)then
          distmin = dist
          intp(1,i) = k
        endif
      enddo
      cntp(1,i) = xtrapfunc(distmin,xtrapdist)
      cntp(2:3,i) = 0d0
    enddo d1
    
    if(numElems > 1000000)then
      write(*,'(2x,A)') '100%'
    endif
    
  endsubroutine

!****************************************************************************
  subroutine interp_coef_tri(xi,yi,xt,yt,w)
! Interpolation coefficients for a linear triangle. 
! Solves the equation for a 2D plane defined by the three points (xt(3),yt(3))
! using Cramer's Rule and calculates the interpolation coefficients w(3) 
! for a point at a point (xi,yi).
!
! Points are normalized from 0 to 1 in order to reduce precision errors.
!
! The method is more expensive than others but has the advantage that the 
! input triangle points (x(3),y(3)) do not need to sorted in any way.
!****************************************************************************
    implicit none
    double precision,intent(in):: xt(3),yt(3),xi,yi
    double precision,intent(out):: w(3)
    double precision :: dInv,xn(3),yn(3),xin,yin,xmin,ymin,xrange,yrange,sumw
    
    !Normalize to reduce precision errors
    xmin = minval(xt)
    ymin = minval(yt)
    xrange = maxval(xt) - xmin
    yrange = maxval(yt) - ymin
    xn = (xt-xmin)/xrange
    yn = (yt-ymin)/yrange
    xin = (xi-xmin)/xrange
    yin = (yi-ymin)/yrange
    
    !Calculate coefficients using Cramer's Rule
    dInv = 1d0/(xn(1)*(yn(2)-yn(3)) + xn(2)*(yn(3)-yn(1)) + xn(3)*(yn(1)-yn(2)))
    w(1) = (xin*(yn(2)-yn(3)) + yin*(xn(3)-xn(2)) + xn(2)*yn(3) - xn(3)*yn(2))*dInv
    w(2) = (xin*(yn(3)-yn(1)) + yin*(xn(1)-xn(3)) + xn(3)*yn(1) - xn(1)*yn(3))*dInv
    w(3) = (xin*(yn(1)-yn(2)) + yin*(xn(2)-xn(1)) + xn(1)*yn(2) - xn(2)*yn(1))*dInv
    
    sumw = sum(w)
    w = w/sum(w)
    if(abs(sumw-1.0)>1.0e-5)then
      write(*,*) 'ERROR: Problem calculating interpolation coefficients'
      write(*,*) 'w = ',w
      write(*,*) 'sum = w(1:3) = ',sumw
      write(*,*) 'xt(1:3) = ',xt
      write(*,*) 'yt(1:3) = ',yt
      write(*,*) 'xi = ',xi
      write(*,*) 'yi = ',yi
      write(*,*) 'd = ',1d0/dInv
      stop
    endif
    
  endsubroutine
  
!************************************************************
  pure function inTriangle(xi,yi,xt,yt,at) result(inTri)
! Determines whether the point (xi,yi) is within the
! triangle (xt(3),yt(3))
! The triangle points can be in any order.
!************************************************************
    implicit none
    double precision,intent(in):: xi,yi,xt(3),yt(3)
    double precision,intent(in),optional:: at
    double precision:: x2(3),y2(3),abc,pab,pbc,pac
    double precision:: suma,err
    logical:: inTri
    
    !Area of ABC triangle
    if(present(at))then
      abc = at
    else  
      abc = triangle_area(xt,yt)
    endif
    
    if(abc <= 1.0e-15)then
      inTri = .false.
      return
    endif
    
    !Initialize
    x2 = xt
    y2 = yt
    
    !Area of PBC triangle
    x2(1) = xi
    y2(1) = yi
    pbc = triangle_area(x2,y2)
    
    !Area of PAC triangle
    x2(2) = xt(1)
    y2(2) = yt(1)
    pac = triangle_area(x2,y2)
    
    !Area of PAB triangle
    x2(3) = xt(2)
    y2(3) = yt(2)
    pab = triangle_area(x2,y2)
    
    !Total area
    suma = pbc + pac + pab
    
    !Normalized error
    err = abs(suma-abc)/abc
    if(err <= 0.0001d0)then
      inTri = .true.
    else
      inTri = .false.
    endif
    
  endfunction

!*************************************************
  pure function triangle_area(xt,yt) result(at)
! Calculates the area of a triangle (x(3),y(3))
! The triangle points can be in any order.
!*************************************************
    implicit none
    double precision,intent(in) :: xt(3),yt(3)
    double precision :: xn(3),yn(3),at
    
    !Subtract min to reduce precision error
    xn = xt - minval(xt)
    yn = yt - minval(yt)
    
    at = 0.5d0*abs(xn(1)*(yn(2)-yn(3)) &
       + xn(2)*(yn(3)-yn(1)) &
       + xn(3)*(yn(1)-yn(2)))
    
  endfunction
  
!*********************************************
  pure function xtrapfunc(dist,xtrapdist)
! Extrapolation function
!*********************************************
    implicit none
    double precision,intent(in) :: dist,xtrapdist
    double precision :: xtrapfunc
    
    xtrapfunc = 0.5 + 0.5*cos(PI*min(dist,xtrapdist)/xtrapdist)
    
  endfunction
  
!****************************************************************************************
  subroutine extract_child_water_levels(wsefilepar,wsefilechl,numNodesChild,nodeParent)
!****************************************************************************************
    implicit none
    character(len=*),intent(in):: wsefilepar,wsefilechl
    integer,intent(in):: numNodesChild !number of child nodes
    integer,intent(in):: nodeParent(numNodesChild)
    double precision:: time, t_inc
    integer:: i,j,numNodesParent,numTimeSteps, step, node, ts_inc, idx
    integer:: numWetNodesParent, numWetNodesChild, ierr, intSteps
    character(len=128) :: input
    integer(kind=8):: intLong
    double precision,allocatable:: elevNodeParent(:)
    double precision:: elev, elevDry
    logical:: isFullFormat = .true.
    logical:: writeStatus = .false.
    double precision,parameter:: ELEV_DRY = -99999.0
    
100 format(e22.10,i15,i10,e17.8)
150 format(f,i,i,f)
200 format(i10,e17.8)
300 format(i11,i11,e16.7,i6,i6)
400 format(a)
    
    !Parent Water Level File
    open(63,file=wsefilepar)
    read(63,'(a80)') input
    read(63,*) numTimeSteps, numNodesParent, t_inc, ts_inc, idx
    allocate(elevNodeParent(numNodesParent))
    
    !Child Water Level File
    open(163,file=wsefilechl)
    write(163,'(a80)') input
    write(163,300) numTimeSteps, numNodesChild, t_inc, ts_inc, idx
    
    intLong = int(numNodesParent,kind=8) * int(numTimeSteps,kind=8)
    if(intLong > 10000000)then
      writeStatus = .true.
      write(*,*) 'Interpolating Subdomain WSE:'
    endif
    intSteps = numTimeSteps/100 + 1
    
    do i=1,numTimeSteps
      !Step Header
      numWetNodesParent = -1
      read(63,150,iostat=ierr) time, step, numWetNodesParent, elevDry !time in seconds
      
      if(numWetNodesParent > 0)then
        isFullFormat = .false.
      else
        elevDry = ELEV_DRY
        isFullFormat = .true.
      endif
      
      if(writeStatus .and. i == 1)then
        if(isFullFormat)then
          write(*,'(2x,A)'), 'Full Format'
        else
          write(*,'(2x,A)'), 'Compact Format'
        endif
      endif
      
      if(writeStatus .and. mod(i,intSteps)==0 .and. i/=numTimeSteps)then
        write(*,'(2x,I2,A)') 100*i/numTimeSteps,'%'
      endif
      
      !Read Parent Node Water Levels
      if(isFullFormat)then !Full format
        do j=1,numNodesParent
          read(63,*) node, elevNodeParent(j)
        enddo
      else !Compact format
        elevNodeParent = ELEV_DRY
        do j=1,numWetNodesParent
          read(63,*) node,elev
          elevNodeParent(node) = elev
        enddo
      endif
      
      !Write Child Node Water Levels
      !if(isFullFormat)then !Full Format
      !  write(163,100) time,step
      !  do j=1,numNodesChild
      !    write(163,200) j, elevNodenumWetNodesChildParent(nodeParent(j))
      !  enddo
      !else !Compact Format
        numWetNodesChild = 0
        do j=1,numNodesChild
          if(abs(elevNodeParent(nodeParent(j)) - ELEV_DRY) > 1e-6)then
            numWetNodesChild = numWetNodesChild + 1
          endif
        enddo
        write(163,100) time, step, numWetNodesChild, ELEV_DRY
        do j=1,numNodesChild
          if(abs(elevNodeParent(nodeParent(j)) - ELEV_DRY) > 1e-6)then
            write(163,200) j, elevNodeParent(nodeParent(j))
          endif
        enddo
      !endif
    enddo
    if(writeStatus)then
      write(*,'(2x,A)') '100%'
    endif
    
    close(63)
    close(163)
    
    !Clean
    deallocate(elevNodeParent)
    
  endsubroutine
  
!*****************************************************
  function inPolygon(np,xp,yp,xi,yi) result(inPoly)
! Determines if the point (xi,yi) is within
! the polygon (xp(np),yp(np)).
! The order of the points does NOT matter.
! If the sum of triangle areas
! is equal to the poly area than the point is
! within the polygon.
!*****************************************************
    implicit none
    integer,intent(in):: np
    double precision,intent(in):: xp(np),yp(np)
    double precision,intent(in):: xi,yi
    double precision:: xpS(np),ypS(np) !Sorted polygon points
    double precision:: ang(np)
    integer:: rank(np)
    logical:: inPoly
    integer:: j,k
    double precision:: xa
    
    !Sort Polygon
    xpS = xp - sum(xp) / dble(np)
    ypS = yp - sum(yp) / dble(np)
    ang = atan2(ypS,xpS)
    call rankArray(ang,rank)
    do k=1,np
      xpS(k) = xp(rank(k))
      ypS(k) = yp(rank(k))
    enddo
    
    inPoly = .false.
    do k=1,np
      if(k < np)then
        j = k + 1
      else
        j = 1
      endif
      xa = (xpS(j) - xpS(k))*(yi - ypS(k)) / (ypS(j) - ypS(k)) + xpS(k)
      if(((((ypS(k) <= yi) .and. (yi < ypS(j))) .or. &
            (ypS(j) <= yi) .and. (yi < ypS(k)))) .and. (xi < xa))then
        inPoly = .not.inPoly
      endif
    enddo
    
  endfunction
  
!**********************************************************************
  pure subroutine rankArray(array,rank)
!Ranks an array in increasing order using the insertion sort algorithm
!The algorithm is very fast compared to other algorithms for 
!arrays smaller than 16 to 20 elements
!**********************************************************************
  implicit none
  double precision,intent(in):: array(:)
  integer,intent(inout):: rank(:) ! inout to avoid automatic deallocation of an allocatable array on entry
  integer:: i,j,tmp,N
  
  N = size(array)
  do i=1,N
    rank(i) = i
  enddo
  do i=2,N
    j = i - 1
    tmp = rank(i)
    do while(j >= 1)
      if(array(rank(j)) < array(tmp))then
        exit
      endif
      rank(j+1) = rank(j)
      j = j - 1
    enddo
    rank(j+1) = tmp
  enddo
  
endsubroutine
  
endmodule