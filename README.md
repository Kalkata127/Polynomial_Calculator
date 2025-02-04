# Polynomial Calculator

This program provides a set of functionalities for working with **polynomials with rational coefficients**.  
It allows users to perform various operations, including arithmetic operations, polynomial evaluation, polynomial division, factorization, and more.  

The calculator implements **mathematical algorithms** such as:
- **Euclidean Algorithm** (for finding the greatest common divisor of polynomials)
- **Horner’s Method** (for efficient polynomial evaluation)
- **Vieta’s Formulas** (for relationships between polynomial roots and coefficients)
- **Binomial Expansion** (for representing polynomials in powers of \(x + a\))
  
### **Available Operations:**
1. **Sum two polynomials**  
2. **Subtract two polynomials**  
3. **Multiply polynomial by a scalar**  
4. **Multiply two polynomials**  
5. **Find the value of a polynomial at a given number**  
6. **Divide two polynomials**  
7. **Compute the GCD of two polynomials**  
8. **Display Vieta’s formulas for a given polynomial**  
9. **Represent a polynomial in powers of \(x + a\)**  
10. **Factor a polynomial and find its rational roots**  

## Usage

### Entering Polynomials  

- **Powers can be entered in any order** and must be separated by spaces. There is no requirement for them to be sequential, meaning some powers can be skipped entirely, as shown in the example below.  
- The program will **always prompt for the coefficient of \( x^0 \)** (the constant term) if it is not included in the initial input.  
- **Coefficients can be rational numbers**, entered in the format `numerator/denominator` (e.g., `3/4`, `-5/2`). If the denominator is `1`, the value can be written as a whole number instead.  
<br>

## Examples: 

### Summing Two Polynomials  
```
Enter your choice: 1
Enter powers for P(x):
Enter powers: 3 1 0
Enter coefficient for power 3: 2
Enter coefficient for power 1: -5/2
Enter coefficient for power 0: 3
Enter powers for Q(x):
Enter powers: 4 3
Enter coefficient for power 4: 1
Enter coefficient for power 3: -2
Enter coefficient for power 0: 7/3
Result of P(x) + Q(x): 1x^4 -5/2x^1 +16/3
```

### Dividing Two Polynomials  
```
Enter your choice: 6
Enter powers for P(x):
Enter powers: 3 2 1
Enter coefficient for power 3: 1
Enter coefficient for power 2: -4/3
Enter coefficient for power 1: -5/6
Enter coefficient for power 0: 1
Enter powers for Q(x):
Enter powers: 2 1
Enter coefficient for power 2: 1
Enter coefficient for power 1: -2/3
Enter coefficient for power 0: -1/2
Quotient: 1x^1 -2/3
Remainder: -7/9x^1 +2/3
```

### GCD of Two Polynomials 
```
Enter your choice: 7
Enter powers for P(x):
Enter powers: 3 2 1
Enter coefficient for power 3: 1
Enter coefficient for power 2: -2
Enter coefficient for power 1: -1
Enter coefficient for power 0: 2
Enter powers for Q(x):
Enter powers: 2 1
Enter coefficient for power 2: 1
Enter coefficient for power 1: 0
Enter coefficient for power 0: -1
GCD of P(x) and Q(x): 1x^2 -1
```

### Vieta's formulas for polynomial
```
Enter your choice: 8
Enter powers: 3 2 1
Enter coefficient for power 3: 2
Enter coefficient for power 2: -5
Enter coefficient for power 1: 4
Enter coefficient for power 0: -3
Calculating Vieta's Formulas for P(x):
P(x) = 2x^3 -5x^2 +4x^1 -3
Vieta's Formulas for polynomial P(x):
x1 + x2 + x3 + ... = 5/2
x1x2 + x1x3 + x2x3 + ... = 2/1
x1x2x3...xn = 3/2
```


### Representing polynomial in powers of (x+a)
```
Enter your choice: 9
Enter powers: 3 2 1
Enter coefficient for power 3: 1
Enter coefficient for power 2: 5/2
Enter coefficient for power 1: 7/4
Enter coefficient for power 0: 9/8
Enter scalar (as a rational number): 3/2
P(x) = 1x^3 +5/2x^2 +7/4x^1 +9/8
Shift value (a): 3/2
P(x + 3/2) = (x + 3/2)^3 + (-2)(x + 3/2)^2 + (x + 3/2) + (3/4)
```
