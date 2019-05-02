def vecScalarMulti(scalar, vector):
  """ Computes scalar vector multiplication 
  Multiplies every element a vector by a scalar and returns the result.

  Args:
    scalar: A number
    vector: A list of numbers representing a vector.
  Returns:
    A list of numbers representing the desired scalar vector multiplication.

  """

  result = [0] * len(vector)
  for iterator in range(len(result)):
    result[iterator] =  scalar * vector[iterator]
  return result  

def degreeFour(x):
  """ Computes the degree 4 vandermonde matrix

        This function takes in a vector, represented as lists, then builds the degree 4 vandermonde matrix A and returns it as a list.
                 
        Args: 
             x: An arbitary vector of arbitary length.
        
        Returns: 
             A degree 4 vandermonde matrix.
  """
  result = []
  for exponent in range(5):
    temp = []
    for element in range(len(x)):
      temp.append(x[element]**exponent)
    result.append(temp)
  return result

def conjugate(z):
  """ Changes the sign of the imaginary part of our number.conjugate

      This function takes a complex number, changes the sign of the imaginary part to get a complex number of different form.

      Args:
            z: A complex number

      Returns:
               A complex number
  """
  result = z.real
  result = result - (z.imag)*1j
  return result

def transpose(A):
  """ Computes the transpose of the matrix

        This function takes in a vandermonde matrix then returns then switches the row and coloumn indices of the matrix to get the transpose.
                 
        Args: 
             A: A vandermonde matrix
        
        Returns: 
             The transpose of the vandermonde matrix A.
  """
  result = []
  for iterator in range(len(A)-1):
    temp = []
    for element in range(len(A)-1):
      temp.append(A[iterator][element])
    result.append(temp)
  return result

def conjugateTranspose(A):
    """ Computes the conjugate transpose of a matrix

        This function takes a matrix in the complex domain and then takes the rows of A and place them as the columns of our outputs and then conjugate the entries.
                 
        Args: 
             A: A matrix in the complex domain
        
        Returns: 
              A matrix in the complex domain.
    """
    result = transpose(A)

    for iterator in range(len(A) -1):
        for element in range(len(A[0]) -1):
          result[iterator][element] = conjugate(result[iterator][element])
    return result

def BackSub(A,b):
  """ Computers the x in the following equation Ax= B

        This function takes a matrix A and a vector b and uses the information for all the following rows to solve for the element of x, corresponding to the current row.
                 
        Args: 
             A: A matrix A and a vector b
        
        Returns: 
              The vector x which solves Ax = B
  """
  for iterator in range(len(A[0])):
    a = len(A[0]-1)
    result = 0
    for i in range(a-iterator+1, len(A)-1):
      result = b(a-iterator - vecScalarMulti(result[i], A[a-iterator][i])* (1/(A[a-iterator]*[a-iterator])))
  return result

def DegreeFourInterpoletion(b,x):
  """ Computes the degree 4 interpoletion

        This function takes in a solution b and returns the degree 4 interpolating polynomial. 
                 
        Args: 
              We have b as the result of the backsubsitution from BackSub function.
        
        Returns: 
              A degree 4 interpolating polynomial.
  """

  return ( b[0] + b[1]**x + b[2]*(x**2) + b[3]*(x**3) + b[4]*(x**4) )
  

def two_norm(vector):
    """ Computes the two_norm of a vector

        This function adds all the elements of the vector and then take the square root of the sum of all the elements of the vector.
                 
        Args: 
              A arbitary vector of arbitary length
        
        Returns: 
              The square root of the sum of all the elements in the vector.
  """
    sum = 0
    y = 0
	
    for i in range(len(vector)):
	    if len(vector) == 0:
			    print("invalid")
	    else:
		    sum = sum + (vector[i])**2
		    y = sum**(0.5)
    return y
		
def dot(vector01, vector02):
    """ Computes the dot product of two vectors

    Takes in two vectors and computes their dot product.

    Args:
      vector01: A list of real numbers representing a vector.
      vector02: A list of real numbers of the same dimensions as vector 01 also reresenting a vector.

    returns: 

      The dot product of the inputs   
  """
    result = 0
    for i in range(len(vector01)):
	    result = result + vector01[i] * vector02[i]
    return result	
	
def ModifiedGramSchmidt(X):
    """ Computes the QR factorization

        This function takes in a matrix and does orthogonal decomposition and normalization.
                 
        Args: 
             We have a unitary matrix Q and upper triangle R.
        
        Returns:
             A m*n matrix A 
  """
    rowsA = len(X)
    colsA = len(X[0])
    Q = [ [0 for row in range(rowsA)] for col in range(colsA) ]
    r = [ [0 for row in range(rowsA)] for col in range(colsA) ]
    for k in range(0,rowsA):
      r[k,k] = two_norm(X[:,k])
      Q[:,k] = (1.0/r[k,k])*X[:,k]
    for j in range(k+1,rowsA):
      r[k,j] = dot(transpose(X[:,j]), Q[:,k])
      X[:,j] = X[:,j] - r[k,j]*Q[:,k]
    return [Q,r]
	