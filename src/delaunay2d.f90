!=================================================================
module delaunay2d
! FORTRAN90 library which computes the Delaunay triangulation 
! of a set of points in the plane. 
! Modified from John Burkardt's TABLE_DELAULAY.f90
!
!  Authors:
!    John Burkardt
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!    Modified by Alex Sanchez for use with HEC-RAS.
!
!  License:
!    This code is distributed under the GNU LGPL license. 
!
!  References:
!    John Burkardt, 
!      TABLE_DELAUNAY Triangulate Points in 2D.
!      https://people.math.sc.edu/Burkardt/f_src/table_delaunay/table_delaunay.html
!    Barry Joe,
!      GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!      Advances in Engineering Software,
!      Volume 13, pages 325-331, 1991.
!=================================================================
implicit none
  
  contains
  
!*****************************************************************************
  subroutine dtris2(point_num, point_xy, tri_num, tri_vert, tri_nabe)
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Input:
!    POINT_NUM : number of vertices
!    point_xy(2,point_num) : coordinates of the vertices
!      On output, the vertices have been sorted into 
!      dictionary order.
!  
!  Output:
!    TRI_NUM : number of triangles in the triangulation
!       tri_num is equal to 2*point_num - NB - 2, where NB is the 
!      number of boundary vertices.
!    TRI_VERT(3,TRI_NUM) :  nodes that make up
!      each triangle.  The elements are indices of POINT_XY.  The vertices of the 
!      triangles are in counter clockwise order.
!    TRI_NABE(3,TRI_NUM) : triangle neighbor list
!      Positive elements are indices of TIL; negative elements are used 
!      for links of a counter clockwise linked list of boundary edges; 
!      LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!      to the neighbor along edge from vertex J to J+1 (mod 3).
!*****************************************************************************
    implicit none
    integer,intent(in) :: point_num
    double precision,intent(inout) :: point_xy(2,point_num)
    integer,intent(out) :: tri_num
    integer,intent(out) :: tri_vert(3,point_num*2)
    integer,intent(out) :: tri_nabe(3,point_num*2)
    double precision :: cmax
    integer :: e, i, ierr, j, k, l, ledg, lr, ltri
    integer :: indx(point_num)
    integer :: m, m1, m2, n, redg, rtri, t, top
    integer :: stack(point_num)
    double precision :: tol
    
    tol = 100D0 * epsilon(tol)
    
    ierr = 0
    
    !Sort the vertices by increasing (x,y).
    call r82vec_sort_heap_index_a(point_num, point_xy, indx)
    
    call r82vec_permute(point_num, indx, point_xy)
    
    !Make sure that the data points are "reasonably" distinct.
    m1 = 1
    
    do i = 2, point_num
      m = m1
      m1 = i
      k = 0
      do j = 1, 2
        cmax = max(abs(point_xy(j,m)), abs(point_xy(j,m1)))
        if(tol * (cmax + 1d0) < abs(point_xy(j,m) - point_xy(j,m1)))then
          k = j
          exit
        endif
      enddo
      if( k == 0 )then
        write(*, '(a)' ) ' '
        write(*, '(a)' ) 'DTRIS2 - Fatal error!'
        write(*, '(a,i8)' ) '  Fails for point number I = ', i
        write(*, '(a,i8)' ) '  M = ', m
        write(*, '(a,i8)' ) '  M1 = ', m1
        write(*, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
        write(*, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
        ierr = 224
        return
      endif
    enddo
    
    !Starting from points M1 and M2, search for a third point M that
    !makes a "healthy" triangle (M1,M2,M)
    m1 = 1
    m2 = 2
    j = 3
    
    do
      if( point_num < j )then
        write(*, '(a)' ) ' '
        write(*, '(a)' ) 'DTRIS2 - Fatal error!'
        ierr = 225
        return
      endif
      m = j
      lr = lrline(point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
        point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0D0)
      if( lr /= 0 )then
        exit
      endif
      j = j + 1
    enddo
    
    !Set up the triangle information for (M1,M2,M), and for any other
    !triangles you created because points were collinear with M1, M2.
    tri_num = j - 2
    
    if( lr == -1 )then
      tri_vert(1,1) = m1
      tri_vert(2,1) = m2
      tri_vert(3,1) = m
      tri_nabe(3,1) = -3
      do i = 2, tri_num
        m1 = m2
        m2 = i+1
        tri_vert(1,i) = m1
        tri_vert(2,i) = m2
        tri_vert(3,i) = m
        tri_nabe(1,i-1) = -3 * i
        tri_nabe(2,i-1) = i
        tri_nabe(3,i) = i - 1
      enddo
      tri_nabe(1,tri_num) = -3 * tri_num - 1
      tri_nabe(2,tri_num) = -5
      ledg = 2
      ltri = tri_num
    else
      tri_vert(1,1) = m2
      tri_vert(2,1) = m1
      tri_vert(3,1) = m
      tri_nabe(1,1) = -4
      do i = 2, tri_num
        m1 = m2
        m2 = i+1
        tri_vert(1,i) = m2
        tri_vert(2,i) = m1
        tri_vert(3,i) = m
        tri_nabe(3,i-1) = i
        tri_nabe(1,i) = -3 * i - 3
        tri_nabe(2,i) = i - 1
      enddo
      tri_nabe(3,tri_num) = -3 * tri_num
      tri_nabe(2,1) = -3 * tri_num - 2
      ledg = 2
      ltri = 1
    endif
    
    !Insert the vertices one at a time from outside the convex hull,
    !determine visible boundary edges, and apply diagonal edge swaps until
    !Delaunay triangulation of vertices (so far) is obtained.
    top = 0
    
    do i = j+1, point_num
      m = i
      m1 = tri_vert(ledg,ltri)
      if( ledg <= 2 )then
        m2 = tri_vert(ledg+1,ltri)
      else
        m2 = tri_vert(1,ltri)
      endif
      lr = lrline(point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
        point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )
      if( 0 < lr )then
        rtri = ltri
        redg = ledg
        ltri = 0
      else
        l = -tri_nabe(ledg,ltri)
        rtri = l / 3
        redg = mod(l,3) + 1
      endif
      
      call vbedg(point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
        tri_vert, tri_nabe, ltri, ledg, rtri, redg)
      
      n = tri_num + 1
      l = -tri_nabe(ledg,ltri)
      
      do
        t = l / 3
        e = mod(l, 3) + 1
        l = -tri_nabe(e,t)
        m2 = tri_vert(e,t)
        if( e <= 2 )then
          m1 = tri_vert(e+1,t)
        else
          m1 = tri_vert(1,t)
        endif
        tri_num = tri_num + 1
        tri_nabe(e,t) = tri_num
        tri_vert(1,tri_num) = m1
        tri_vert(2,tri_num) = m2
        tri_vert(3,tri_num) = m
        tri_nabe(1,tri_num) = t
        tri_nabe(2,tri_num) = tri_num - 1
        tri_nabe(3,tri_num) = tri_num + 1
        top = top + 1
        if( point_num < top )then
          ierr = 8
          write(*, '(a)' ) ' '
          write(*, '(a)' ) 'DTRIS2 - Fatal error!'
          write(*, '(a)' ) '  Stack overflow.'
          return
        endif
        stack(top) = tri_num
        if( t == rtri .and. e == redg )then
          exit
        endif
      enddo
      
      tri_nabe(ledg,ltri) = -3 * n - 1
      tri_nabe(2,n) = -3 * tri_num - 2
      tri_nabe(3,tri_num) = -l
      ltri = n
      ledg = 2
      
      call swapec(m, top, ltri, ledg, point_num, point_xy, tri_num, &
        tri_vert, tri_nabe, stack, ierr)
      
      if( ierr /= 0 )then
        write(*, '(a)' ) ' '
        write(*, '(a)' ) 'DTRIS2 - Fatal error!'
        write(*, '(a)' ) '  Error return from SWAPEC.'
        return
      endif
    enddo
    
    !Now account for the sorting that we did.
    do i = 1, 3
      do j = 1, tri_num
        tri_vert(i,j) = indx(tri_vert(i,j))
      enddo
    enddo
    
    call perm_inverse(point_num, indx)
    
    call r82vec_permute(point_num, indx, point_xy)
    
  endsubroutine
  
!*****************************************************************************
  function diaedg(x0, y0, x1, y1, x2, y2, x3, y3)
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Input:
!    x0, y0, x1, y1, x2, y2, x3, y3 : coordinates of the vertices of 
!      a quadrilateral, given in counter clockwise order.
!  
!  Output:
!    DIAEDG : chooses a diagonal:
!      +1, if diagonal edge 02 is chosen;
!      -1, if diagonal edge 13 is chosen;
!       0, if the four vertices are cocircular.
!*****************************************************************************
    implicit none
    double precision,intent(in) :: x0, x1, x2, x3, y0, y1, y2, y3
    integer :: diaedg
    double precision :: ca, cb
    double precision :: dx10, dx12, dx30, dx32, dy10, dy12, dy30, dy32
    double precision :: s, tol, tola, tolb
    
    tol = 100.0D+00 * epsilon(tol)
    
    dx10 = x1 - x0
    dy10 = y1 - y0
    dx12 = x1 - x2
    dy12 = y1 - y2
    dx30 = x3 - x0
    dy30 = y3 - y0
    dx32 = x3 - x2
    dy32 = y3 - y2

    tola = tol * max(abs(dx10), abs(dy10), abs(dx30), abs(dy30))
    tolb = tol * max(abs(dx12), abs(dy12), abs(dx32), abs(dy32))

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if( tola < ca .and. tolb < cb )then
      diaedg = -1
    elseif( ca < -tola .and. cb < -tolb )then
      diaedg = 1
    else
      tola = max(tola, tolb)
      s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca
      if( tola < s )then
        diaedg = -1
      elseif( s < -tola )then
        diaedg = 1
      else
        diaedg = 0
      endif
    endif
    
  endfunction

!*****************************************************************************
  function i4_modp(i, j) result(value)
!! i4_modp returns the nonnegative remainder of I4 division.
!
!  Discussion:
!    If
!      nRem = i4_modp(i, j)
!      nMult = ( I - nRem ) / j
!   then
!      i = j * nMult + nRem
!    where nRem is always nonnegative.
!
!    The mod function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, i4_modp(A,360) is between 0 and 360, always.
!
!  Example:
!        i     j     MOD I4_MODP    Factorization
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Input:
!    i : number to be divided
!    j : number that divides i
!  
!  Output:
!    I4_MODP :nonnegative remainder when i is divided by j.
!*****************************************************************************
    implicit none
    integer,intent(in) :: i, j
    integer :: value
    
    if( j == 0 )then
      write(*, '(a)' ) ' '
      write(*, '(a)' ) 'I4_MODP - Fatal error!'
      write(*, '(a,i8)' ) '  Illegal divisor J = ', j
      stop
    endif
    
    value = mod(i, j)
    
    if( value < 0 )then
      value = value + abs(j)
    endif
    
  end function

!*****************************************************************************
  function i4_sign(x)
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!    An I4 is an integer :: value.
!
!  Input:
!    X : number whose sign is desired
!
!  Output:
!    I4_SIGN : sign of X
!*****************************************************************************
    implicit none
    integer,intent(in) :: x
    integer :: i4_sign
    
    if( x < 0 )then
      i4_sign = -1
    else
      i4_sign = +1
    endif
    
  endfunction

!*****************************************************************************
  function i4_wrap(ival, ilo, ihi)
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!    An I4 is an integer :: value.
!
!  Example:
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Input:
!    IVAL : value
!    ILO, IHI : desired bounds
!  
!  Output:
!   I4_WRAP : a "wrapped" version of the value
!*****************************************************************************
    implicit none
    integer,intent(in) :: ival, ilo, ihi
    integer :: i4_wrap, value, wide, jlo, jhi
    
    jlo = min(ilo, ihi)
    jhi = max(ilo, ihi)
    
    wide = jhi - jlo + 1
    
    if( wide == 1 )then
      value = jlo
    else
      value = jlo + i4_modp(ival - jlo, wide)
    endif
    
    i4_wrap = value
    
  endfunction

!*****************************************************************************
  function lrline(xu, yu, xv1, yv1, xv2, yv2, dv)
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Input:
!    XU, YU : coordinates of the point whose
!      position relative to the directed line is to be determined.
!    XV1, YV1, XV2, YV2 :  coordinates of two points
!      that determine the directed base line.
!    DV : signed distance of the directed line
!      from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!      DV is positive for a line to the left of the base line.
!
!  Output:
!    LRLINE : result:
!      +1, the point is to the right of the directed line;
!       0, the point is on the directed line;
!      -1, the point is to the left of the directed line.
!*****************************************************************************
    implicit none
    double precision,intent(in) :: xu, yu, xv1, yv1, xv2, yv2, dv
    double precision :: dx, dxu, dy, dyu, t, tol, tolabs
    integer :: lrline
    
    tol = 100.0D+00 * epsilon(tol)
    
    dx = xv2 - xv1
    dy = yv2 - yv1
    dxu = xu - xv1
    dyu = yu - yv1
    
    tolabs = tol * max(abs(dx), abs(dy), abs(dxu), abs(dyu), abs(dv))
    
    t = dy * dxu - dx * dyu + dv * sqrt(dx * dx + dy * dy)
    
    if( tolabs < t )then
      lrline = +1 !Right
    elseif( -tolabs <= t )then
      lrline = 0 !Directed
    else
      lrline = -1 !Left
    endif
    
  endfunction

!*****************************************************************************
  subroutine perm_check(n, p, base, ierror)
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Input:
!    n: number of entries
!    p(N) : array to check
!    base : index base
!
!  Output:
!    ierror : error flag
!      0, the array represents a permutation.
!      nonzero, the array does not represent a permutation.
!      The smallest missing value is equal to IERROR.
!*****************************************************************************
    implicit none
    integer,intent(in) :: n, base
    integer,intent(in) :: p(n)
    integer,intent(out) :: ierror
    integer :: f, seek
    
    ierror = 0
    
    do seek = base, base + n - 1
      ierror = 1
      do f = 1, n
        if( p(f) == seek )then
          ierror = 0
          exit
        endif
      enddo
      if( ierror /= 0 )then
        write(*, '(a)' ) ' '
        write(*, '(a)' ) 'PERM_CHECK - Fatal error!'
        write(*, '(a)' ) '  The input array does not represent'
        write(*, '(a)' ) '  a proper permutation.'
        stop
      endif
    enddo
    
  endsubroutine

!*****************************************************************************
  subroutine perm_inverse(n, p)
!! PERM_INVERSE inverts a permutation "in place".
!
!  Input:
!    N :  number of objects being permuted.
!  Input/Output:
!    P(N) :permutation, in standard index form.
!      On output, P describes the inverse permutation
!*****************************************************************************
    implicit none
    integer,intent(in) :: n
    integer,intent(inout) :: p(n)
    integer, parameter :: base = 1
    integer :: i, i0, i1, i2, ierror, is
    
    if( n <= 0 )then
      write(*, '(a)' ) ' '
      write(*, '(a)' ) 'PERM_INVERSE - Fatal error!'
      write(*, '(a,i8)' ) '  Input value of N = ', n
      stop
    endif
    
    call perm_check(n, p, base, ierror)
    
    if( ierror /= 0 )then
      write(*, '(a)' ) ' '
      write(*, '(a)' ) 'PERM_INVERSE - Fatal error!'
      write(*, '(a)' ) '  PERM_CHECK rejects this permutation.'
      stop
    endif
    
    is = 1
    
    do i = 1, n
      i1 = p(i)
      do while ( i < i1 )
        i2 = p(i1)
        p(i1) = -i2
        i1 = i2
      enddo
      is = - i4_sign ( p(i) )
      p(i) = is * abs ( p(i) )
    enddo
    
    do i = 1, n
      i1 = - p(i)
      if( 0 <= i1 )then
        i0 = i
        do
          i2 = p(i1)
          p(i1) = i0
          if( i2 < 0 )then
            exit
          endif
          i0 = i1
          i1 = i2
        enddo
      endif
    enddo
    
  endsubroutine

!*****************************************************************************
  subroutine r82vec_permute(n, p, a)
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!    An R82VEC is an array of pairs of R8 values.
!
!    The same logic can be used to permute an array of objects of any 
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!    Input:
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Input:
!    N : number of objects
!    P(N) : permutation array
!      P(I) = J means that the I-th element of the output array should 
!      be the J-th element of the input array
!  
!  Input/output:
!    A(2,N) : array to be permuted
!*****************************************************************************
    implicit none
    integer,intent(in) :: n
    integer,intent(inout) :: p(n)
    integer, parameter :: dim_num = 2
    double precision,intent(inout) :: a(dim_num,n)
    double precision :: a_temp(dim_num)
    integer, parameter :: base = 1
    integer :: ierror, iget, iput, istart
    
    call perm_check(n, p, base, ierror)
    
    if( ierror /= 0 )then
      write(*, '(a)' ) ' '
      write(*, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
      write(*, '(a)' ) '  PERM_CHECK rejects this permutation.'
      stop
    endif
    
    !Search for the next element of the permutation that has not been used.
    do istart = 1, n
      if( p(istart) < 0 )then
        cycle
      elseif( p(istart) == istart )then
        p(istart) = - p(istart)
        cycle
      else
        a_temp(1:dim_num) = a(1:dim_num,istart)
        iget = istart

        !Copy the new value into the vacated entry.
        do
          iput = iget
          iget = p(iget)
          p(iput) = - p(iput)
          if( iget < 1 .or. n < iget )then
            write(*, '(a)' ) ' '
            write(*, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
            write(*, '(a)' ) '  A permutation index is out of range.'
            write(*, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
            stop
          endif
          if( iget == istart )then
            a(1:dim_num,iput) = a_temp(1:dim_num)
            exit
          endif
          a(1:dim_num,iput) = a(1:dim_num,iget)
        enddo
      endif
    enddo
    
    !Restore the signs of the entries.
    p(1:n) = - p(1:n)
    
  end subroutine

!*****************************************************************************
  subroutine r82vec_sort_heap_index_a(n, a, indx)
!! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!
!  Discussion:
!    An R82VEC is an array of R82's.
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly":
!      A(1:2,INDX(1:N)) is sorted,
!    or explicitly, by the call
!      call r82vec_permute( n, indx, a)
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Input:
!    n : number of entries in the array
!    a(2,n) : array to be index-sorted
! 
!  Output:
!    indx(n) : sort index
!      The I-th element of the sorted array is A(1:2,INDX(I))
!*****************************************************************************
    implicit none
    integer,intent(in) :: n
    integer, parameter :: dim_num = 2
    double precision, intent(in) :: a(dim_num,n)
    integer, intent(out) :: indx(n)
    integer :: i, indxt, ir, j, l
    double precision :: aval(dim_num)
    
    if( n < 1 )then
      return
    endif
    
    do i = 1, n
      indx(i) = i
    enddo
    
    if( n == 1 )then
      return
    endif
    
    l = n / 2 + 1
    ir = n
    
    do
      if( 1 < l )then
        l = l - 1
        indxt = indx(l)
        aval(1:dim_num) = a(1:dim_num,indxt)
      else
        indxt = indx(ir)
        aval(1:dim_num) = a(1:dim_num,indxt)
        indx(ir) = indx(1)
        ir = ir - 1
        if( ir == 1 )then
          indx(1) = indxt
          exit
        endif
      endif
      
      i = l
      j = l + l
      
      do while ( j <= ir )
        if( j < ir )then
          if(   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
               ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
                 a(2,indx(j)) <  a(2,indx(j+1)) ) )then
            j = j + 1
          endif
        endif
        if(   aval(1) <  a(1,indx(j)) .or. &
             ( aval(1) == a(1,indx(j)) .and. &
               aval(2) <  a(2,indx(j)) ) )then
          indx(i) = indx(j)
          i = j
          j = j + j
        else
          j = ir + 1
        endif
      enddo
      indx(i) = indxt
    enddo
    
  endsubroutine

!*****************************************************************************
  subroutine swapec(i, top, btri, bedg, point_num, point_xy, tri_num, &
    tri_vert, tri_nabe, stack, ierr)
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Input:
!    I : index of the new vertex
!    POINT_NUM, the number of points
!    POINT_XY(2,POINT_NUM) : the coordinates of the points
!    TRI_NUM : number of triangles
!
!  Input/Output:
!    TOP, index of the top of the stack
!      On output, TOP is zero
!    BTRI, BEDG; on input, if positive, are 
!      the triangle and edge indices of a boundary edge whose updated indices
!      must be recorded.  On output, these may be updated because of swaps.
!    TRI_VERT(3,TRI_NUM) : triangle incidence list
!      May be updated on output because of swaps.
!    TRI_NABE(3,TRI_NUM) : triangle eighbor list
!      negative values are used for links of the counter-clockwise 
!      linked list of boundary edges;  May be updated on output because of swaps.
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!    STACK(MAXST) : Workspace. on input, entries 1 through
!      TOP contain the indices of initial triangles (involving vertex I)
!      put in stack; the edges opposite I should be in interior;  entries
!      TOP+1 through MAXST are used as a stack.
!  
!  Output:
!    IERR : set to 8 for abnormal return.
!*****************************************************************************
    implicit none
    integer,intent(in):: i
    integer,intent(inout):: top, bedg, btri
    integer,intent(in) :: point_num, tri_num
    double precision,intent(in) :: point_xy(2,point_num)
    integer,intent(inout) :: tri_nabe(3,tri_num)
    integer,intent(inout) :: tri_vert(3,tri_num)
    integer,intent(inout) :: stack(point_num)
    integer :: a, b, c, e, ee, em1, ep1
    integer :: f, fm1, fp1, ierr, l, r, s
    integer :: swap, t, tt, u
    double precision :: x, y
    
    !Determine whether triangles in stack are Delaunay, and swap
    !diagonal edge of convex quadrilateral if not.
    x = point_xy(1,i)
    y = point_xy(2,i)
    
    do
      if( top <= 0 )then
        exit
      endif
      t = stack(top)
      top = top - 1
      if( tri_vert(1,t) == i )then
        e = 2
        b = tri_vert(3,t)
      elseif( tri_vert(2,t) == i )then
        e = 3
        b = tri_vert(1,t)
      else
        e = 1
        b = tri_vert(2,t)
      endif
      a = tri_vert(e,t)
      u = tri_nabe(e,t)
      if( tri_nabe(1,u) == t )then
        f = 1
        c = tri_vert(3,u)
      elseif( tri_nabe(2,u) == t )then
        f = 2
        c = tri_vert(1,u)
      else
        f = 3
        c = tri_vert(2,u)
      endif
      
      swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
        point_xy(2,c), point_xy(1,b), point_xy(2,b) )
      
      if( swap == 1 )then
        em1 = i4_wrap ( e - 1, 1, 3 )
        ep1 = i4_wrap ( e + 1, 1, 3 )
        fm1 = i4_wrap ( f - 1, 1, 3 )
        fp1 = i4_wrap ( f + 1, 1, 3 )
        
        tri_vert(ep1,t) = c
        tri_vert(fp1,u) = i
        r = tri_nabe(ep1,t)
        s = tri_nabe(fp1,u)
        tri_nabe(ep1,t) = u
        tri_nabe(fp1,u) = t
        tri_nabe(e,t) = s
        tri_nabe(f,u) = r
        
        if( 0 < tri_nabe(fm1,u) )then
          top = top + 1
          stack(top) = u
        endif
        if( 0 < s )then
          if( tri_nabe(1,s) == u )then
            tri_nabe(1,s) = t
          elseif( tri_nabe(2,s) == u )then
            tri_nabe(2,s) = t
          else
            tri_nabe(3,s) = t
          endif
          top = top + 1
          if( point_num < top )then
            ierr = 8
            return
          endif
          stack(top) = t
        else
          if( u == btri .and. fp1 == bedg )then
            btri = t
            bedg = e
          endif
          l = - ( 3 * t + e - 1 )
          tt = t
          ee = em1
          do while ( 0 < tri_nabe(ee,tt) )
            tt = tri_nabe(ee,tt)
            if( tri_vert(1,tt) == a )then
              ee = 3
            elseif( tri_vert(2,tt) == a )then
              ee = 1
            else
              ee = 2
            endif
          enddo
          tri_nabe(ee,tt) = l
        endif
        
        if( 0 < r )then
          if( tri_nabe(1,r) == t )then
            tri_nabe(1,r) = u
          elseif( tri_nabe(2,r) == t )then
            tri_nabe(2,r) = u
          else
            tri_nabe(3,r) = u
          endif
        else
          if( t == btri .and. ep1 == bedg )then
            btri = u
            bedg = f
          endif
          l = - ( 3 * u + f - 1 )
          tt = u
          ee = fm1
          do while ( 0 < tri_nabe(ee,tt) )
            tt = tri_nabe(ee,tt)
            if( tri_vert(1,tt) == b )then
              ee = 3
            elseif( tri_vert(2,tt) == b )then
              ee = 1
            else
              ee = 2
            endif
          enddo
          tri_nabe(ee,tt) = l
        endif
      endif
    enddo
    
  endsubroutine
  
!*****************************************************************************
  subroutine vbedg(x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
    ltri, ledg, rtri, redg)
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Input:
!    x, y : coordinates of a point outside
!      the convex hull of the current triangulation
!    point_num : number of points
!    point_xy : coordinates of the vertices
!    tri_num : number of triangles
!    tri_vert :  triangle incidence list
!    tri_nabe : triangle neighbor list; negative values are used for links of 
!      a counter clockwise linked list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!   
!  Input/output:
!    ltri, ledg.  If ltri /= 0 then these 
!      values are assumed to be already computed and are not changed, else they 
!      are updated.  On output, LTRI is the index of boundary triangle to the 
!      left of the leftmost boundary triangle visible from (X,Y), and LEDG is 
!      the boundary edge of triangle LTRI to the left of the leftmost boundary
!      edge visible from (X,Y).  1 <= LEDG <= 3.
!    rtri : On input, the index of the 
!      boundary triangle to begin the search at.  On output, the index of the 
!      rightmost boundary triangle visible from (X,Y).
!    redg : edge of triangle RTRI that 
!      is visible from (X,Y).  1 <= REDG <= 3.
!*****************************************************************************
    implicit none
    integer,intent(in) :: point_num, tri_num
    integer,intent(inout) :: ledg, ltri, redg, rtri
    double precision,intent(in) :: point_xy(2,point_num)
    integer,intent(in) :: tri_nabe(3,tri_num)
    integer,intent(in) :: tri_vert(3,tri_num)
    double precision :: x, y
    integer :: a, b, e, l, lr, t
    logical :: ldone
    
    !Find the rightmost visible boundary edge using links,then possibly
    !leftmost visible boundary edge using triangle neighbor information.
    if( ltri == 0 )then
      ldone = .false.
      ltri = rtri
      ledg = redg
    else
      ldone = .true.
    endif
    
    do
      l = -tri_nabe(redg,rtri)
      t = l / 3
      e = mod(l, 3) + 1
      a = tri_vert(e,t)
      if( e <= 2 )then
        b = tri_vert(e+1,t)
      else
        b = tri_vert(1,t)
      endif
      lr = lrline(x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), point_xy(2,b), 0D0)
      if( lr <= 0 )then
        exit
      endif
      rtri = t
      redg = e
    enddo
    
    if( ldone )then
      return
    endif
    
    t = ltri
    e = ledg
    
    do
      b = tri_vert(e,t)
      e = i4_wrap ( e-1, 1, 3 )
      do while ( 0 < tri_nabe(e,t) )
        t = tri_nabe(e,t)
        if( tri_vert(1,t) == b )then
          e = 3
        elseif( tri_vert(2,t) == b )then
          e = 1
        else
          e = 2
        endif
      enddo
      
      a = tri_vert(e,t)
      
      lr = lrline(x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), point_xy(2,b), 0D0)
      
      if( lr <= 0 )then
        exit
      endif
    enddo
    
    ltri = t
    ledg = e
    
  endsubroutine
  
endmodule