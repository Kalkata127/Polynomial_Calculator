/**
*
* Solution to course project # 4
* Introduction to programming course
* Faculty of Mathematics and Informatics of Sofia University
* Winter semester 2023/2024
*
* @author Kalin Simeonov
* @idnumber 8MI0600459* @compiler VC
*
* <Main cpp file>
*
*/
#include <iostream>
#include <vector>

using namespace std;

const int MAX_ARR_SIZE = 50;

// ---------------------- Utils -------------------------------------
struct Rational {
    int numerator;
    int denominator;
};

Rational createRational(int numerator, int denominator = 1) {
    Rational r;
    r.numerator = numerator;
    r.denominator = denominator;
    return r;
}

void sortVector(vector<pair<int, Rational>>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        for (size_t j = 0; j < vec.size() - i - 1; j++) {
            if (vec[j].first < vec[j + 1].first) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

int GCD(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return abs(a);
}

int LCM(int a, int b) {
    if (a == 0 || b == 0) return 0;
    return abs(a * b) / GCD(a, b);
}

Rational BinomialCoefficient(int n, int k) {
    if (k < 0 || k > n) return createRational(0, 1); 
    if (k == 0 || k == n) return createRational(1, 1); 

    Rational result = createRational(1, 1);

    for (int i = 1; i <= k; i++) {
        result = multiplyRationals(result, createRational(n - i + 1, 1));
        result = simplifyFraction(result); 

        result = divideRationals(result, createRational(i, 1));
        result = simplifyFraction(result);
    }

    return result;
}


// ---------------------- Rational Operations -----------------------
Rational addRationals(const Rational& r1, const Rational& r2) {
    int numerator = r1.numerator * r2.denominator + r2.numerator * r1.denominator;
    int denominator = r1.denominator * r2.denominator;
    int gcd = GCD(numerator, denominator);
    return createRational(numerator / gcd, denominator / gcd);
}

Rational subtractRationals(const Rational& r1, const Rational& r2) {
    int numerator = r1.numerator * r2.denominator - r2.numerator * r1.denominator;
    int denominator = r1.denominator * r2.denominator;
    int gcd = GCD(numerator, denominator);
    return createRational(numerator / gcd, denominator / gcd);
}

Rational multiplyRationals(const Rational& r1, const Rational& r2) {
    int numerator = r1.numerator * r2.numerator;
    int denominator = r1.denominator * r2.denominator;
    int gcd = GCD(numerator, denominator);
    return createRational(numerator / gcd, denominator / gcd);
}

Rational divideRationals(const Rational& r1, const Rational& r2) {
    if (r2.numerator == 0) {
        cout << "Error: Division by zero. Returning default rational (0/1)." << endl;
        return createRational(0, 1);
    }

    int numerator = r1.numerator * r2.denominator;
    int denominator = r1.denominator * r2.numerator;

    if (denominator == 0) {
        cout << "Error: Invalid rational number. Denominator cannot be zero. Returning default rational (0/1)." << endl;
        return createRational(0, 1);
    }

    int gcd = GCD(abs(numerator), abs(denominator));
    numerator /= gcd;
    denominator /= gcd;

    return createRational(numerator, denominator);
}

Rational sqrtRational(const Rational& r) {
    if (r.numerator < 0 || r.denominator <= 0) {
        cout << "Error: Square root of a negative or invalid rational number is undefined. Returning default rational (0/1)." << endl;
        return createRational(0, 1);
    }

    int sqrtNumerator = sqrt(r.numerator);
    int sqrtDenominator = sqrt(r.denominator);

    if (sqrtNumerator * sqrtNumerator == r.numerator && sqrtDenominator * sqrtDenominator == r.denominator) {
        return createRational(sqrtNumerator, sqrtDenominator);
    }
    else {
        cout << "Error: Result of square root is not a rational number. Returning default rational (0/1)." << endl;
        return createRational(0, 1);
    }
}

Rational powRational(const Rational& base, int exponent) {
    if (exponent == 0) return createRational(1, 1);

    Rational result = base;
    for (int i = 1; i < exponent; i++) {
        result = multiplyRationals(result, base);
    }
    return result;
}

Rational AddRational(const Rational& r1, const Rational& r2) {
    int lcm = LCM(r1.denominator, r2.denominator);
    int numerator1 = r1.numerator * (lcm / r1.denominator);
    int numerator2 = r2.numerator * (lcm / r2.denominator);

    int sumNumerator = numerator1 + numerator2;
    int gcd = GCD(sumNumerator, lcm);

    return createRational(sumNumerator / gcd, lcm / gcd);
}

Rational simplifyFraction(const Rational& r) {
   int gcd = GCD(abs(r.numerator), abs(r.denominator));
   return { r.numerator / gcd, r.denominator / gcd };
}


// ---------------------- Displayment -------------------------------
void DisplayPolynomial(const vector<pair<int, Rational>> &vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        const pair<int, Rational> term = vec[i];
        int power = term.first;
        Rational coeff = term.second;

        if (coeff.numerator == 0) continue;

        if (coeff.numerator > 0 && i > 0) {
            cout << "+";
        }
        cout << coeff.numerator;
        if (coeff.denominator != 1) cout << "/" << coeff.denominator;

        if (power != 0) cout << "x^" << power << " ";
    }
    cout << endl;
}

void DisplayRoots(const vector<pair<Rational, int>>& factors) {
    cout << "RATIONAL ROOTS:" << endl;
    for (size_t i = 0; i < factors.size(); i++) {
        PrintRoot(factors[i].first);
        cout << " -> " << factors[i].second << "-fold root" << endl;
    }
}

void PrintRational(const Rational& r) {
    cout << r.numerator;
    if (r.denominator != 1) {
        cout << "/" << r.denominator;
    }
}

void PrintPolynomial(const vector<pair<int, Rational>>& polynomial, const Rational& shift) {
    cout << "P(x + "; PrintRational(shift); cout << ") = ";

    bool firstTerm = true;
    for (size_t i = 0; i < polynomial.size(); i++) {
        if (polynomial[i].second.numerator == 0) continue;

        if (!firstTerm) cout << " + ";
        firstTerm = false;

        if (polynomial[i].second.numerator != 1 || polynomial[i].second.denominator != 1 || polynomial[i].first == 0) {
            cout << "("; PrintRational(polynomial[i].second); cout << ")";
        }

        if (polynomial[i].first > 0) {
            cout << "(x + "; PrintRational(shift); cout << ")";
            if (polynomial[i].first > 1) cout << "^" << polynomial[i].first;
        }
    }
    cout << endl;
}

void displayFactors(const vector<pair<Rational, int>>& factors) {
    for (size_t i = 0; i < factors.size(); i++) {
        const pair<Rational, int>& factor = factors[i];
        Rational root = factor.first;
        int multiplicity = factor.second;

        if (i > 0) {
            cout << " * ";
        }

        cout << "(x ";
        if (root.numerator < 0) {
            cout << "+ " << -root.numerator;
        }
        else {
            cout << "- " << root.numerator;
        }

        if (root.denominator != 1 && root.denominator != -1) {
            cout << "/" << abs(root.denominator);
        }
        cout << ")";

        if (multiplicity > 1) {
            cout << "^" << multiplicity;
        }
    }
    cout << endl;
}

void PrintRoot(const Rational& root) {
    cout << "x = ";
    if (root.denominator == -1) cout << -root.numerator;
    else if (root.denominator == 1) cout << root.numerator;
    else cout << root.numerator << "/" << root.denominator;
}

// ---------------------- Input Handling ----------------------------
void getInput(char* input) {
    if (input == nullptr) return;
    bool validInput = false;

    while (!validInput) {
        int index = 0;
        char c;

        while (index < MAX_ARR_SIZE - 1) {
            c = cin.get();

            if (c == '\n') {
                break;
            }

            input[index++] = c;
        }

        input[index] = '\0';

        if (index == MAX_ARR_SIZE - 1 && c != '\n') {
            while (cin.get() != '\n');
            cout << "Input was too long. Please try again." << endl;
        }
        else {
            validInput = true;
        }
    }
}

void EnterCoefficients(const int* powers, int powerCount, vector<pair<int, Rational>>& coefficients, char* input) {
    if (powers == nullptr || input == nullptr) {
        cout << "Error: Null pointer passed to EnterCoefficients." << endl;
        return;
    }

    coefficients.clear();

    for (int i = 0; i < powerCount; i++) {
        bool valid = false;
        Rational coefficient;

        do {
            cout << "Enter coefficient for power " << powers[i] << ": ";
            getInput(input);

            if (isCoefficientValid(input)) {
                ProcessRational(input, coefficient);
                valid = true;
            }
            else {
                cout << "Invalid coefficient, please re-enter." << endl;
            }
        } while (!valid);

        if (coefficient.numerator != 0) {
            coefficients.emplace_back(powers[i], coefficient);
        }
    }
}

void EnterPowers(int* powers, char* input, int& counter) {
    if (powers == nullptr || input == nullptr) {
        cout << "Error: Null pointer passed to EnterPowers." << endl;
        return;
    }

    do {
        cout << "Enter powers: ";
        getInput(input);
        counter = 0;
    } while (!arePowersValid(input));

    ProcessPowers(input, powers, counter);
    counter = RemoveDuplicates(powers, counter);

    bool zeroPowerFound = false;
    for (int i = 0; i < counter; ++i) {
        if (powers[i] == 0) {
            zeroPowerFound = true;
            break;
        }
    }
    if (!zeroPowerFound) {
        powers[counter] = 0; 
        ++counter;            
    }
}

Rational InputScalar(char* input) {
    Rational scalar;

    while (true) {
        cout << "Enter scalar (as a rational number): ";
        getInput(input);
        if (isCoefficientValid(input) && ProcessRational(input, scalar)) {
            break;
        }
        cout << "Invalid scalar. Please re-enter." << endl;
    }
    return scalar;
}

void InputPolynomial(vector<pair<int, Rational>>& poly, int* powers, char* input, int& counter) {
    EnterPowers(powers, input, counter);
    EnterCoefficients(powers, counter, poly, input);
    sortVector(poly);
}

int splitInput(const char* input, int* output) {
    if (!input || !output) return 0;

    int index = 0, num = 0;
    bool isNeg = false, inNum = false;

    while (*input != '\0') {
        if (*input == '-') {
            if (inNum) return -1;
            isNeg = true;
        }
        else if (isDigit(*input)) {
            num = num * 10 + (*input - '0');
            inNum = true;
        }
        else if (*input == ' ' && inNum) {
            output[index++] = isNeg ? -num : num;
            num = 0, isNeg = false, inNum = false;
        }
        input++;
    }

    if (inNum) {
        output[index++] = isNeg ? -num : num;
    }

    return index;
}


// ---------------------- Input Validating --------------------------
bool isDigit(char c) {
    return c >= '0' && c <= '9';
}

bool isInputValid(const char input[]) {
    if (input == nullptr) {
        cout << "Error: Null pointer passed as input!" << endl;
        return false;
    }

    if (input[0] == '\0') {
        cout << "Input is empty!" << endl;
        return false;
    }

    bool hasNonSpaceCharacter = false;

    for (int i = 0; input[i] != '\0'; i++) {
        char c = input[i];

        if (!(c == '-' || c == '/' || isDigit(c) || c == ' ')) {
            cout << "Input contains forbidden characters!" << endl;
            return false;
        }

        if (c != ' ') {
            hasNonSpaceCharacter = true;
        }
    }

    if (!hasNonSpaceCharacter) {
        cout << "Input can't be empty or contain only spaces!" << endl;
        return false;
    }

    return true;
}

bool arePowersValid(const char input[]) {
    const int Int_max = 2147483647; 
    long long power = 0;            
    bool inNumber = false;

    for (int i = 0; input[i] != '\0'; i++) {
        char c = input[i];

        if (isDigit(c)) {            
            power = power * 10 + (c - '0');
            if (power > Int_max) {
                cout << "Input contains a power that is too large!" << endl;
                return false;
            }
            inNumber = true;
        }
        else if (c == ' ' && inNumber) {     
            power = 0;
            inNumber = false;
        }
        else if (c != ' ') {                 
            cout << "Input contains invalid power format!" << endl;
            return false;
        }
    }

    return true;
}

bool isCoefficientValid(const char input[]) {
    if (input[0] == '0' && input[1] == '\0') return true;
    if (!isInputValid(input)) return false;

    bool hasSlash = false, hasMinus = false;
    int i = 0, denominator = 0;

    if (input[i] == '-') {
        hasMinus = true;
        i++;
    }

    bool inNumber = false;
    for (i; input[i] != '\0'; i++) {
        char c = input[i];

        if (c == '/') {
            if (hasSlash || !inNumber) return false;
            hasSlash = true;
            inNumber = false;
        }
        else if (isDigit(c)) {
            inNumber = true;
            if (hasSlash) denominator = denominator * 10 + (c - '0');
        }
        else return false;
    }

    return inNumber && (!hasSlash || denominator != 0);
}

bool getValidatedNumber(const char* input, int& number) {
    for (size_t i = 0; input[i] != '\0'; i++) {
        if (i >= 2 || input[i] < '0' || input[i] > '9') { 
            cout << "Invalid input."<<endl;
            return false;
        }
    }

    number = 0;
    for (size_t i = 0; input[i] != '\0'; i++) {
        number = number * 10 + (input[i] - '0');
    }

    return true;
}


// ---------------------- Processing --------------------------------
void ProcessPowers(char* input, int* pPowers, int& counter) {
    if (input == nullptr || pPowers == nullptr) {
        return;
    }
    counter = splitInput(input, pPowers);

    if (counter == -1) {
        cout << "Error processing input for powers" << endl;
        return;
    }
}

bool ProcessRational(const char* input, Rational& rational) {
    int numerator = 0, denominator = 0;
    bool isNegative = false, inDenominator = false;

    if (input[0] == '0' && input[1] == '\0') {
        rational = createRational(0, 1);
        return true; 
    }

    for (int i = 0; input[i] != '\0'; i++) {
        char c = input[i];

        if (c == '-') {
            if (i == 0 || (inDenominator && denominator == 0)) isNegative = true;
            else return false;
        }
        else if (c == '/') {
            if (inDenominator || numerator == 0) return false;
            inDenominator = true;
        }
        else if (isDigit(c)) {
            if (inDenominator) {
                denominator = denominator * 10 + (c - '0');
            }
            else {
                numerator = numerator * 10 + (c - '0');
            }
        }
        else return false;
    }

    if (inDenominator && denominator == 0) return false;

    if (numerator == 0) {
        rational = createRational(0, 1);
        return true;
    }

    if (isNegative) numerator = -numerator;
    rational = createRational(numerator, inDenominator ? denominator : 1);
    return true;
}

int RemoveDuplicates(int* powers, int count) {
    int uniqueCount = 0;

    for (int i = 0; i < count; i++) {
        bool isDuplicate = false;

        for (int j = 0; j < uniqueCount; j++) {
            if (powers[i] == powers[j]) {
                isDuplicate = true;
                break;
            }
        }

        if (!isDuplicate) {
            powers[uniqueCount++] = powers[i];
        }
    }

    return uniqueCount;
}


// ---------------------- Operations --------------------------------
void PolynomSum(vector<pair<int, Rational>>* v1Ptr, vector<pair<int, Rational>>* v2Ptr, vector<pair<int, Rational>>& resVec) {
    if (!v1Ptr || !v2Ptr) {
        cout << "Error: Null pointer passed to PolynomSum." << endl;
        return;
    }

    const vector<pair<int, Rational>>& vec1 = *v1Ptr;
    const vector<pair<int, Rational>>& vec2 = *v2Ptr;

    size_t idV1 = 0, idV2 = 0;

    while (idV1 < vec1.size() || idV2 < vec2.size()) {
        if (idV1 < vec1.size() && (idV2 >= vec2.size() || vec1[idV1].first > vec2[idV2].first)) {
            resVec.push_back(vec1[idV1++]);
        }
        else if (idV2 < vec2.size() && (idV1 >= vec1.size() || vec2[idV2].first > vec1[idV1].first)) {
            resVec.push_back(vec2[idV2++]);
        }
        else {
            Rational sumCoeff = AddRational(vec1[idV1].second, vec2[idV2].second);
            resVec.emplace_back(vec1[idV1].first, sumCoeff);
            idV1++;
            idV2++;
        }
    }
}

void PolynomSubtract(vector<pair<int, Rational>>* v1Ptr, vector<pair<int, Rational>>* v2Ptr, vector<pair<int, Rational>>& resVec) {
    if (!v1Ptr || !v2Ptr) {
        cout << "Error: Null pointer passed to PolynomSubtract." << endl;
        return;
    }

    const vector<pair<int, Rational>>& vec1 = *v1Ptr;
    const vector<pair<int, Rational>>& vec2 = *v2Ptr;

    size_t idV1 = 0, idV2 = 0;
    resVec.clear();

    while (idV1 < vec1.size() || idV2 < vec2.size()) {
        if (idV1 < vec1.size() && (idV2 >= vec2.size() || vec1[idV1].first > vec2[idV2].first)) {
            resVec.push_back(vec1[idV1++]);
        }
        else if (idV2 < vec2.size() && (idV1 >= vec1.size() || vec2[idV2].first > vec1[idV1].first)) {
            resVec.emplace_back(vec2[idV2].first, createRational(-vec2[idV2].second.numerator, vec2[idV2].second.denominator));
            idV2++;
        }
        else {
            Rational diffCoeff = AddRational(vec1[idV1].second, createRational(-vec2[idV2].second.numerator, vec2[idV2].second.denominator));
            if (diffCoeff.numerator != 0) {
                resVec.emplace_back(vec1[idV1].first, diffCoeff);
            }
            idV1++;
            idV2++;
        }
    }
}

void HandlePolynomialAddition(vector<pair<int, Rational>>& pVector, vector<pair<int, Rational>>& qVector, vector<pair<int, Rational>>& result, int* pPowers, int* qPowers, char* input, int& pCounter, int& qCounter) {
    cout << "Enter powers for P(x):" << endl;
    InputPolynomial(pVector, pPowers, input, pCounter);

    cout << "Enter powers for Q(x):" << endl;
    InputPolynomial(qVector, qPowers, input, qCounter);

    PolynomSum(&pVector, &qVector, result);
    cout << "Result of P(x) + Q(x): ";
    DisplayPolynomial(result);
}

void HandlePolynomialSubtraction(vector<pair<int, Rational>>& pVector, vector<pair<int, Rational>>& qVector, vector<pair<int, Rational>>& result, int* pPowers, int* qPowers, char* input, int& pCounter, int& qCounter) {
    cout << "Enter powers for P(x):" << endl;
    InputPolynomial(pVector, pPowers, input, pCounter);

    cout << "Enter powers for Q(x):" << endl;
    InputPolynomial(qVector, qPowers, input, qCounter);

    PolynomSubtract(&pVector, &qVector, result);
    cout << "Result of P(x) - Q(x): ";
    DisplayPolynomial(result);
}

void MultiplyByScalar(const Rational& scalar, vector<pair<int, Rational>>& polynomial) {
    if (scalar.numerator == 0) {
        for (size_t i = 0; i < polynomial.size(); i++) {
            polynomial[i].second.numerator = 0; 
            polynomial[i].second.denominator = 1; 
        }
        return;
    }

    for (size_t i = 0; i < polynomial.size(); i++) {
        polynomial[i].second.numerator *= scalar.numerator;

        polynomial[i].second.denominator *= scalar.denominator;

        int greatestCommonDivisor = GCD(abs(polynomial[i].second.numerator), abs(polynomial[i].second.denominator));
        polynomial[i].second.numerator /= greatestCommonDivisor;
        polynomial[i].second.denominator /= greatestCommonDivisor;

        if (polynomial[i].second.denominator < 0) {
            polynomial[i].second.numerator = -polynomial[i].second.numerator;
            polynomial[i].second.denominator = -polynomial[i].second.denominator;
        }
    }
}

void MultiplyPolynomials(const vector<pair<int, Rational>>& p1, const vector<pair<int, Rational>>& p2, vector<pair<int, Rational>>& result) {
    result.clear();
    if (p1.empty() || p2.empty()) return;

    vector<pair<int, Rational>> tempResult;

    for (size_t i = 0; i < p1.size(); i++) {
        for (size_t j = 0; j < p2.size(); j++) {
            int power = p1[i].first + p2[j].first;
            Rational coefficient = createRational(
                p1[i].second.numerator * p2[j].second.numerator,
                p1[i].second.denominator * p2[j].second.denominator
            );

            if (coefficient.denominator < 0) {
                coefficient.numerator = -coefficient.numerator;
                coefficient.denominator = -coefficient.denominator;
            }

            bool combined = false;
            for (size_t k = 0; k < tempResult.size(); k++) {
                if (tempResult[k].first == power) {
                    tempResult[k].second = AddRational(tempResult[k].second, coefficient);
                    combined = true;
                    break;
                }
            }

            if (!combined) {
                tempResult.emplace_back(power, coefficient);
            }
        }
    }

    for (size_t i = 0; i < tempResult.size(); i++) {
        if (tempResult[i].second.numerator != 0) {
            result.emplace_back(tempResult[i]);
        }
    }

    sortVector(result);
}

void EvaluatePolynomial(const vector<pair<int, Rational>>& polynomial, const Rational& input, Rational& result) {
    result = createRational(0, 1);
    Rational powerValue = createRational(1, 1);

    for (size_t i = 0; i < polynomial.size(); i++) {
        int power = polynomial[i].first;
        const Rational& coefficient = polynomial[i].second;

        Rational xPower = createRational(1, 1);
        for (int j = 0; j < power; j++) {
            xPower.numerator *= input.numerator;
            xPower.denominator *= input.denominator;
        }

        Rational term;
        term.numerator = coefficient.numerator * xPower.numerator;
        term.denominator = coefficient.denominator * xPower.denominator;

        result.numerator = result.numerator * term.denominator + term.numerator * result.denominator;
        result.denominator *= term.denominator;

        int gcd = GCD(abs(result.numerator), abs(result.denominator));
        result.numerator /= gcd;
        result.denominator /= gcd;

        if (result.denominator < 0) {
            result.numerator = -result.numerator;
            result.denominator = -result.denominator;
        }
    }
}

void DividePolynomials(const vector<pair<int, Rational>>& P1, const vector<pair<int, Rational>>& P2, vector<pair<int, Rational>>& Q, vector<pair<int, Rational>>& R) {
    Q.clear();
    R = P1;

    while (!R.empty() && R.front().first >= P2.front().first) {
        int powerDiff = R.front().first - P2.front().first;
        Rational coeffQuotient = divideRationals(R.front().second, P2.front().second);

        Q.emplace_back(powerDiff, coeffQuotient);

        vector<pair<int, Rational>> tempPoly;
        for (size_t i = 0; i < P2.size(); i++) {
            tempPoly.emplace_back(
                P2[i].first + powerDiff,
                createRational(P2[i].second.numerator * coeffQuotient.numerator, P2[i].second.denominator * coeffQuotient.denominator)
            );
        }

        vector<pair<int, Rational>> newR;
        size_t i = 0, j = 0;
        while (i < R.size() || j < tempPoly.size()) {
            if (i < R.size() && (j >= tempPoly.size() || R[i].first > tempPoly[j].first)) {
                newR.emplace_back(R[i].first, R[i].second);
                i++;
            }
            else if (j < tempPoly.size() && (i >= R.size() || tempPoly[j].first > R[i].first)) {
                newR.emplace_back(tempPoly[j].first, createRational(-tempPoly[j].second.numerator, tempPoly[j].second.denominator));
                j++;
            }
            else {
                Rational diff;
                diff.numerator = R[i].second.numerator * tempPoly[j].second.denominator - tempPoly[j].second.numerator * R[i].second.denominator;
                diff.denominator = R[i].second.denominator * tempPoly[j].second.denominator;

                if (diff.numerator != 0) {
                    int gcd = GCD(abs(diff.numerator), abs(diff.denominator));
                    diff.numerator /= gcd;
                    diff.denominator /= gcd;
                    newR.emplace_back(R[i].first, diff);
                }
                i++;
                j++;
            }
        }
        R = newR;
    }
}

void PolynomialGCD(const vector<pair<int, Rational>>& P, const vector<pair<int, Rational>>& Q, vector<pair<int, Rational>>& GCD) {
    vector<pair<int, Rational>> A = P;
    vector<pair<int, Rational>> B = Q;
    vector<pair<int, Rational>> quotient, remainder;

    while (!B.empty()) {
        DividePolynomials(A, B, quotient, remainder);
        A = B;
        B = remainder;
    }

    if (!A.empty()) {
        Rational leadingCoeff = A.front().second;
        for (auto& term : A) {
            term.second = divideRationals(term.second, leadingCoeff);
        }
    }

    GCD = A;
}

void VietaFormulas(const vector<pair<int, Rational>>& polynomial) {
    if (polynomial.empty()) {
        cout << "The polynomial is empty!" << endl;
        return;
    }

    cout << "P(x) = ";
    DisplayPolynomial(const_cast<vector<pair<int, Rational>>&>(polynomial));

    cout << "Vieta's Formulas for polynomial P(x):" << endl;

    vector<pair<int, Rational>> sortedPolynomial = polynomial;
    sortVector(sortedPolynomial);

    int degree = sortedPolynomial[0].first;    
    Rational leadingCoeff = sortedPolynomial[0].second; 

    for (size_t i = 1; i < sortedPolynomial.size(); i++) {
        int powerDiff = degree - sortedPolynomial[i].first;
        Rational coeff = sortedPolynomial[i].second;

        Rational term = divideRationals(
            createRational(pow(-1, powerDiff) * coeff.numerator, coeff.denominator),
            leadingCoeff
        );

        if (powerDiff == 1) {
            cout << "x1 + x2 + x3 + ... = " << term.numerator << "/" << term.denominator << endl;
        }
        else if (powerDiff == 2) {
            cout << "x1x2 + x1x3 + x2x3 + ... = " << term.numerator << "/" << term.denominator << endl;
        }
        else if (powerDiff == sortedPolynomial.size() - 1) {
            cout << "x1x2x3...xn = " << term.numerator << "/" << term.denominator << endl;
        }
        else {
            cout << "Sum of products of roots taken " << powerDiff << " at a time: "
                << term.numerator << "/" << term.denominator << endl;
        }
    }
}

void RepresentPolynomialInPowersOfXPlusA(vector<pair<int, Rational>>& polynomial, const Rational& shift) {
    if (polynomial.empty()) {
        cout << "The polynomial is empty!" << endl;
        return;
    }

    cout << "P(x) = "; DisplayPolynomial(polynomial);
    cout << "Shift value (a): "; PrintRational(shift); cout << endl;

    vector<pair<int, Rational>> expandedPolynomial;
    for (size_t i = 0; i < polynomial.size(); i++) {
        ExpandPolynomial(expandedPolynomial, polynomial[i].first, polynomial[i].second, shift);
    }

    sortVector(expandedPolynomial);
    PrintPolynomial(expandedPolynomial, shift);
}

void FactorizePolynomial(vector<pair<int, Rational>>& polynomial) {
    if (polynomial.empty()) {
        cout << "The polynomial is empty!" << endl;
        return;
    }

    cout << "Factoring Polynomial P(x): ";
    DisplayPolynomial(polynomial);

    vector<pair<Rational, int>> factors;
    vector<pair<int, Rational>> remainder = polynomial;

    ExtractFactors(remainder, factors);

    if (remainder.size() == 2) {
        handleQuadratic(remainder);
    }

    DisplayRoots(factors);

    cout << "Factored Polynomial: ";
    displayFactors(factors);
}


// ---------------------- Factorization & Root Finding --------------
void ExtractFactors(vector<pair<int, Rational>>& remainder, vector<pair<Rational, int>>& factors) {
    while (remainder.size() > 1) {
        Rational root;
        if (!findRoot(remainder, root)) {
            cout << "Remaining polynomial is irreducible over the rationals: ";
            DisplayPolynomial(remainder);
            return;
        }

        int multiplicity = 0;
        vector<pair<int, Rational>> divisor = {
            {1, createRational(1, 1)},
            {0, multiplyRationals(createRational(-1, 1), root)}
        };

        vector<pair<int, Rational>> quotient, newRemainder;

        do {
            DividePolynomials(remainder, divisor, quotient, newRemainder);
            Rational valueAtRoot;
            EvaluatePolynomial(remainder, root, valueAtRoot);

            if (valueAtRoot.numerator == 0) {
                ++multiplicity;
                remainder = quotient;
            }
            else {
                break;
            }
        } while (true);

        factors.emplace_back(root, multiplicity);
    }
}

void ExpandPolynomial(vector<pair<int, Rational>>& expandedPolynomial, int power, Rational coeff, const Rational& shift) {
    for (int k = 0; k <= power; k++) {
        Rational binomialCoeff = BinomialCoefficient(power, k);
        Rational shiftTerm = powRational(shift, power - k);
        if ((power - k) % 2 != 0) shiftTerm.numerator = -shiftTerm.numerator;

        Rational termCoeff = multiplyRationals(multiplyRationals(binomialCoeff, shiftTerm), coeff);

        bool termFound = false;
        for (size_t j = 0; j < expandedPolynomial.size(); j++) {
            if (expandedPolynomial[j].first == k) {
                expandedPolynomial[j].second = addRationals(expandedPolynomial[j].second, termCoeff);
                termFound = true;
                break;
            }
        }
        if (!termFound) expandedPolynomial.emplace_back(k, termCoeff);
    }
}

vector<int> findDivisors(int num) {
    vector<int> divisors;

    if (num == 0) {
        cout << "Warning: Zero has no divisors." << endl;
        return divisors;
    }

    num = abs(num);

    for (int i = 1; i <= num; i++) {
        if (num % i == 0) {
            divisors.emplace_back(i);
            divisors.emplace_back(-i);
        }
    }

    return divisors;
}

vector<Rational> generatePossibleRoots(const vector<pair<int, Rational>>& polynomial) {
    if (polynomial.empty()) return {};

    int leadingCoeff = polynomial[0].second.numerator;
    int constantTerm = polynomial.back().second.numerator;

    vector<int> divisorsConstant = findDivisors(constantTerm);
    vector<int> divisorsLeading = findDivisors(leadingCoeff);

    vector<Rational> roots;
    for (size_t i = 0; i < divisorsConstant.size(); i++) {
        for (size_t j = 0; j < divisorsLeading.size(); j++) {
            if (divisorsLeading[j] != 0) {
                Rational fraction = { divisorsConstant[i], divisorsLeading[j] };
                roots.push_back(simplifyFraction(fraction));
            }
        }
    }

    vector<Rational> uniqueRoots;
    for (size_t i = 0; i < roots.size(); i++) {
        bool isDuplicate = false;
        for (size_t j = 0; j < uniqueRoots.size(); j++) {
            if (roots[i].numerator == uniqueRoots[j].numerator &&
                roots[i].denominator == uniqueRoots[j].denominator) {
                isDuplicate = true;
                break;
            }
        }
        if (!isDuplicate) {
            uniqueRoots.push_back(roots[i]);
        }
    }

    return uniqueRoots;
}

bool findRoot(const vector<pair<int, Rational>>& polynomial, Rational& root) {
    vector<Rational> possibleRoots = generatePossibleRoots(polynomial);

    for (size_t i = 0; i < possibleRoots.size(); i++) {
        Rational valueAtCandidate;
        EvaluatePolynomial(polynomial, possibleRoots[i], valueAtCandidate);

        if (valueAtCandidate.numerator == 0) {
            root = possibleRoots[i];
            return true;
        }
    }
    return false;
}

void handleQuadratic(const vector<pair<int, Rational>>& quadratic) {
    if (quadratic.size() != 2) {
        cout << "Error: Expected a quadratic polynomial." << endl;
        return;
    }

    Rational a = quadratic[0].second;
    Rational b = quadratic[1].second;
    Rational c = quadratic[2].second;

    Rational discriminant = subtractRationals(
        multiplyRationals(b, b),
        multiplyRationals(createRational(4, 1), multiplyRationals(a, c))
    );

    if (discriminant.numerator < 0) {
        cout << "Quadratic polynomial has complex roots: ";
        DisplayPolynomial(quadratic);
        return;
    }

    Rational sqrtDiscriminant = sqrtRational(discriminant);

    Rational root1 = divideRationals(
        addRationals(multiplyRationals(createRational(-1, 1), b), sqrtDiscriminant),
        multiplyRationals(createRational(2, 1), a)
    );

    Rational root2 = divideRationals(
        subtractRationals(multiplyRationals(createRational(-1, 1), b), sqrtDiscriminant),
        multiplyRationals(createRational(2, 1), a)
    );

    cout << "Quadratic roots: x = " << root1.numerator << "/" << root1.denominator;
    cout << " and x = " << root2.numerator << "/" << root2.denominator << endl;
}


// ---------------------- Main ---------------------------------------
void ClearAll(vector<pair<int, Rational>>& pVector, vector<pair<int, Rational>>& qVector, vector<pair<int, Rational>>& result, vector<pair<int, Rational>>& quotient, vector<pair<int, Rational>>& remainder, int& pCounter, int& qCounter, int* pPowers, int* qPowers) {
    pVector.clear();
    qVector.clear();
    result.clear();
    quotient.clear();
    remainder.clear();
    pCounter = 0;
    qCounter = 0;

    fill(pPowers, pPowers + MAX_ARR_SIZE, 0);
    fill(qPowers, qPowers + MAX_ARR_SIZE, 0);
}

void DisplayOptions() {
    cout << "\nOptions:" << endl;
    cout << "1 - Sum two polynomials" << endl;
    cout << "2 - Subtract two polynomials" << endl;
    cout << "3 - Multiply polynomial by scalar" << endl;
    cout << "4 - Multiply two polynomials" << endl;
    cout << "5 - Find value of polynomial at a given number" << endl;
    cout << "6 - Divide two polynomials" << endl;
    cout << "7 - GCD of Polynomials" << endl;
    cout << "8 - Display Vieta's formulas for a given polynomial" << endl;
    cout << "9 - Represent a polynomial in powers of (x+a)" << endl;
    cout << "10 - Factor polynomial and find its rational roots" << endl;
    cout << "11 - Quit the application" << endl;
}

int main() {
    cout << "Welcome to the Polynomial Calculator!" << endl;

    int pPowers[MAX_ARR_SIZE] = { 0 }, qPowers[MAX_ARR_SIZE] = { 0 };
    char input[MAX_ARR_SIZE];
    int pCounter = 0, qCounter = 0;

    vector<pair<int, Rational>> pVector, qVector;
    vector<pair<int, Rational>> result;
    vector<pair<int, Rational>> quotient, remainder;

    while (true) {
        DisplayOptions();

        int choice = 0;
        bool validChoice = false;

        do {
            cout << "Enter your choice: ";
            getInput(input);
            validChoice = getValidatedNumber(input, choice);
        } while (!validChoice);

        switch (choice) {
        case 1:
            HandlePolynomialAddition(pVector, qVector, result, pPowers, qPowers, input, pCounter, qCounter);
            break;

        case 2:
            HandlePolynomialSubtraction(pVector, qVector, result, pPowers, qPowers, input, pCounter, qCounter);
            break;

        case 3: {
            Rational scalar = InputScalar(input);
            InputPolynomial(pVector, pPowers, input, pCounter);

            MultiplyByScalar(scalar, pVector);
            cout << "Result of P(x) multiplied by " << scalar.numerator << "/" << scalar.denominator << ": ";
            DisplayPolynomial(pVector);
            break;
        }

        case 4: {
            cout << "Enter powers for P(x):" << endl;
            InputPolynomial(pVector, pPowers, input, pCounter);

            cout << "Enter powers for Q(x):" << endl;
            InputPolynomial(qVector, qPowers, input, qCounter);

            MultiplyPolynomials(pVector, qVector, result);
            cout << "Result of P(x) * Q(x): ";
            DisplayPolynomial(result);
            break;
        }

        case 5: {
            InputPolynomial(pVector, pPowers, input, pCounter);
            Rational inputValue = InputScalar(input);

            result.clear();
            Rational resultValue;
            EvaluatePolynomial(pVector, inputValue, resultValue);

            result.emplace_back(0, resultValue);
            cout << "Result of evaluating P(x) at x = " << inputValue.numerator << "/" << inputValue.denominator << ": ";
            DisplayPolynomial(result);
            break;
        }

        case 6: {
            cout << "Enter powers for P(x):" << endl;
            InputPolynomial(pVector, pPowers, input, pCounter);

            cout << "Enter powers for Q(x):" << endl;
            InputPolynomial(qVector, qPowers, input, qCounter);

            DividePolynomials(pVector, qVector, quotient, remainder);

            cout << "Quotient: ";
            DisplayPolynomial(quotient);

            cout << "Remainder: ";
            DisplayPolynomial(remainder);
            break;
        }

        case 7: {
            cout << "Enter powers for P(x):" << endl;
            InputPolynomial(pVector, pPowers, input, pCounter);

            cout << "Enter powers for Q(x):" << endl;
            InputPolynomial(qVector, qPowers, input, qCounter);

            PolynomialGCD(pVector, qVector, result);
            cout << "GCD of P(x) and Q(x): ";
            DisplayPolynomial(result);
            break;
        }

        case 8:
            InputPolynomial(pVector, pPowers, input, pCounter);
            cout << "Calculating Vieta's Formulas for P(x):" << endl;
            VietaFormulas(pVector);
            break;

        case 9:
            InputPolynomial(pVector, pPowers, input, pCounter);
            Rational shiftValue = InputScalar(input);
            RepresentPolynomialInPowersOfXPlusA(pVector, shiftValue);
            break;

        case 10:
            InputPolynomial(pVector, pPowers, input, pCounter);
            FactorizePolynomial(pVector);
            break;

        case 11:
            cout << "Thank you for using the Polynomial Calculator. Goodbye!" << endl;
            return 0;

        default:
            cout << "Invalid choice. Please enter available ones:" << endl;
            break;
        }

        ClearAll(pVector, qVector, result, quotient, remainder, pCounter, qCounter, pPowers, qPowers);
    }

    return 0;
}
