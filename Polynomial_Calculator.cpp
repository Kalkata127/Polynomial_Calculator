#include <iostream>
#include <vector>

using namespace std;

const int MAX_ARR_SIZE = 50;

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

void DisplayPolynomial(const vector<pair<int, Rational>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        const auto& term = vec[i];
        int power = term.first;
        const Rational& coeff = term.second;

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
        throw invalid_argument("Division by zero: second rational number's numerator is zero.");
    }

    int numerator = r1.numerator * r2.denominator;
    int denominator = r1.denominator * r2.numerator;

    if (denominator == 0) {
        throw invalid_argument("Invalid rational number: denominator cannot be zero.");
    }

    int gcd = GCD(abs(numerator), abs(denominator));
    numerator /= gcd;
    denominator /= gcd;

    return createRational(numerator, denominator);
}

Rational sqrtRational(const Rational& r) {
    if (r.numerator < 0 || r.denominator <= 0) {
        throw invalid_argument("Square root of a negative or invalid rational number is undefined.");
    }

    int sqrtNumerator = sqrt(r.numerator);
    int sqrtDenominator = sqrt(r.denominator);

    if (sqrtNumerator * sqrtNumerator == r.numerator && sqrtDenominator * sqrtDenominator == r.denominator) {
        return createRational(sqrtNumerator, sqrtDenominator);
    }
    else {
        throw invalid_argument("Result of square root is not rational.");
    }
}


bool isInputValid(const char input[]) {
    if ((input == nullptr) ? (cout << "Error: Null pointer passed as input!" << endl, true) : false) return false;

    if ((input[0] == '\0') ? (cout << "Input is empty!" << endl, true) : false) return false;

    bool hasNonSpaceCharacter = false;

    for (int i = 0; input[i] != '\0'; i++) {
        char c = input[i];

        if (!(c == '-' || c == '/' || (c >= '0' && c <= '9') || c == ' ')) {
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

        if (c >= '0' && c <= '9') {            
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
    if (input[0] == '0' && input[1] == '\0') {
        return true;
    }
    if (!isInputValid(input)) {
        return false;
    }

    bool hasSlash = false, hasMinus = false, inNumber = false;
    int length = 0, denominator = 0;

    for (int i = 0; input[i] != '\0'; i++, length++) {
        char c = input[i];

        if (c == '-') {
            if (i != 0 || hasMinus) {
                cout << "Invalid coefficient: '-' can only appear at the beginning!" << endl;
                return false;
            }
            hasMinus = true;
        }
        else if (c == '/') {
            if (hasSlash) {
                cout << "Invalid coefficient: More than one '/' is not allowed!" << endl;
                return false;
            }
            if (!inNumber) {
                cout << "Invalid coefficient: No numerator before '/'!" << endl;
                return false;
            }
            hasSlash = true;
            inNumber = false;
        }
        else if (c >= '0' && c <= '9') {
            if (hasSlash) {
                denominator = denominator * 10 + (c - '0');
            }
            inNumber = true;
        }
        else {
            cout << "Invalid character in coefficient!" << endl;
            return false;
        }
    }

    if (length > 0 && (input[length - 1] == '/' || input[length - 1] == '-')) {
        cout << "Invalid coefficient: Cannot end with '/' or '-'" << endl;
        return false;
    }

    if (hasSlash && denominator == 0) {
        cout << "Invalid coefficient: Denominator cannot be zero!" << endl;
        return false;
    }

    return true;
}


int splitInput(const char* input, int* output) {
    if (input == nullptr || output == nullptr) {
        return 0;
    }

    int index = 0, num = 0;
    bool isNeg = false, inNum = false;

    while (*input != '\0') {
        if (*input == '-') {
            if (inNum) {
                cout << "Invalid input: multiple '-' in a single number" << endl;
                return -1;
            }
            isNeg = true;
        }
        else if (*input >= '0' && *input <= '9') {
            num = num * 10 + (*input - '0');
            inNum = true;
        }
        else if (*input == ' ') {
            if (inNum) {
                output[index++] = isNeg ? -num : num;
                num = 0;
                isNeg = false;
                inNum = false;
            }
        }
        input++;
    }

    if (inNum) {
        output[index++] = isNeg ? -num : num;
    }

    return index;
}

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

int RemoveDuplicates(int* powers, int count) {
    int uniqueCount = 0;

    for (int i = 0; i < count; ++i) {
        bool isDuplicate = false;

        for (int j = 0; j < uniqueCount; ++j) {
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

bool ProcessRational(const char* input, Rational& rational) {
    int numerator = 0;
    int denominator = 0;
    bool isNegative = false;
    bool inDenominator = false;
    bool hasNumerator = false;
    bool hasDenominator = false;

    for (int i = 0; input[i] != '\0'; i++) {
        if (input[i] == '-') {
            if (i == 0 || (inDenominator && denominator == 0)) {
                isNegative = true;
            }
            else {
                return false;
            }
        }
        else if (input[i] == '/') {
            if (inDenominator) {
                return false;
            }
            inDenominator = true;
        }
        else if (input[i] >= '0' && input[i] <= '9') {
            if (inDenominator) {
                denominator = denominator * 10 + (input[i] - '0');
                hasDenominator = true;
            }
            else {
                numerator = numerator * 10 + (input[i] - '0');
                hasNumerator = true;
            }
        }
        else {
            return false;
        }
    }

    if (!hasNumerator) {
        return false; 
    }

    if (inDenominator && !hasDenominator) {
        return false;
    }

    if (inDenominator && denominator == 0) {
        return false;
    }

    if (isNegative) {
        numerator = -numerator;
    }

    if (numerator == 0) {
        rational = createRational(0, 1);
        return true;
    }

    rational = createRational(numerator, inDenominator ? denominator : 1);
    return true;
}


void EnterCoefficients(const int* powers, int powerCount, vector<pair<int, Rational>>& coefficients, char* input) {
    if (powers == nullptr || input == nullptr) {
        cout << "Error: Null pointer passed to EnterCoefficients." << endl;
        return;
    }

    for (int i = 0; i < powerCount; ++i) {
        bool valid = false;

        do {
            cout << "Enter coefficient for power " << powers[i] << ": ";
            getInput(input);

            if (isCoefficientValid(input)) {
                Rational coefficient;
                ProcessRational(input, coefficient);
                coefficients.emplace_back(powers[i], coefficient);
                valid = true;
            }
            else {
                cout << "Invalid coefficient. Please re-enter." << endl; 
            }
        } while (!valid);
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

Rational AddRational(const Rational& r1, const Rational& r2) {
    int lcm = LCM(r1.denominator, r2.denominator);
    int numerator1 = r1.numerator * (lcm / r1.denominator);
    int numerator2 = r2.numerator * (lcm / r2.denominator);

    int sumNumerator = numerator1 + numerator2;
    int gcd = GCD(sumNumerator, lcm);

    return createRational(sumNumerator / gcd, lcm / gcd);
}

void PolynomSum(vector<pair<int, Rational>>* v1Ptr, vector<pair<int, Rational>>* v2Ptr, vector<pair<int, Rational>>& resVec) {
    if (v1Ptr == nullptr || v2Ptr == nullptr) {
        cout << "Error: Null pointer passed to PolynomSum." << endl;
        return;
    }

    vector<pair<int, Rational>>& vec1 = *v1Ptr, & vec2 = *v2Ptr;

    size_t idV1 = 0, idV2 = 0;

    while (idV1 < vec1.size() && idV2 < vec2.size()) {
        pair<int, Rational> term1 = vec1[idV1], term2 = vec2[idV2];

        if (term1.first == term2.first) { 
            Rational sumCoeff = AddRational(term1.second, term2.second);
            resVec.emplace_back(term1.first, sumCoeff);
            idV1++;
            idV2++;
        }
        else if (term1.first > term2.first) { 
            resVec.emplace_back(term1.first, term1.second);
            idV1++;
        }
        else {
            resVec.emplace_back(term2.first, term2.second);
            idV2++;
        }
    }

    while (idV1 < vec1.size()) {
        pair<int, Rational> remainingTerm = vec1[idV1]; 
        resVec.emplace_back(remainingTerm.first, remainingTerm.second);
        idV1++;
    }

    while (idV2 < vec2.size()) {
        pair<int, Rational> remainingTerm = vec2[idV2]; 
        resVec.emplace_back(remainingTerm.first, remainingTerm.second);
        idV2++;
    }
}

void PolynomSubtract(vector<pair<int, Rational>>* v1Ptr, vector<pair<int, Rational>>* v2Ptr, vector<pair<int, Rational>>& resVec) {
    if (v1Ptr == nullptr || v2Ptr == nullptr) {
        cout << "Error: Null pointer passed to PolynomSubtract." << endl;
        return;
    }

    vector<pair<int, Rational>>& vec1 = *v1Ptr;
    vector<pair<int, Rational>>& vec2 = *v2Ptr;

    size_t idV1 = 0, idV2 = 0;

    resVec.clear();

    while (idV1 < vec1.size() && idV2 < vec2.size()) {
        pair<int, Rational> term1 = vec1[idV1], term2 = vec2[idV2];

        if (term1.first == term2.first) {
            Rational diffCoeff = AddRational(term1.second, createRational(-term2.second.numerator, term2.second.denominator));
            if (diffCoeff.numerator != 0) {
                resVec.emplace_back(term1.first, diffCoeff);
            }
            idV1++;
            idV2++;
        }
        else if (term1.first > term2.first) {
            resVec.emplace_back(term1.first, term1.second);
            idV1++;
        }
        else {
            Rational negCoeff = createRational(-term2.second.numerator, term2.second.denominator);
            resVec.emplace_back(term2.first, negCoeff);
            idV2++;
        }
    }

    while (idV1 < vec1.size()) {
        resVec.emplace_back(vec1[idV1].first, vec1[idV1].second);
        idV1++;
    }

    while (idV2 < vec2.size()) {
        Rational negCoeff = createRational(-vec2[idV2].second.numerator, vec2[idV2].second.denominator);
        resVec.emplace_back(vec2[idV2].first, negCoeff);
        idV2++;
    }
}

void MultiplyByScalar(const Rational& scalar, vector<pair<int, Rational>>& polynomial) {
    if (scalar.numerator == 0) {
        for (size_t i = 0; i < polynomial.size(); ++i) {
            polynomial[i].second.numerator = 0; 
            polynomial[i].second.denominator = 1; 
        }
        return;
    }

    for (size_t i = 0; i < polynomial.size(); ++i) {
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

void MultiplyPolynomials(const vector<pair<int, Rational>>& p1,const vector<pair<int, Rational>>& p2,vector<pair<int, Rational>>& result)
{
    result.clear();

    if (p1.empty() || p2.empty()) {
        return;
    }

    vector<pair<int, Rational>> tempResult;

    for (size_t i = 0; i < p1.size(); i++) {
        for (size_t j = 0; j < p2.size(); j++) {
            Rational coefficient = createRational(
                p1[i].second.numerator * p2[j].second.numerator,
                p1[i].second.denominator * p2[j].second.denominator
            );


            if (coefficient.denominator < 0) {
                coefficient.numerator = -coefficient.numerator;
                coefficient.denominator = -coefficient.denominator;
            }

            int power = p1[i].first + p2[j].first;

            tempResult.push_back({ power, coefficient });
        }
    }

    for (size_t i = 0; i < tempResult.size(); i++) {
        bool combined = false;
        for (size_t j = 0; j < result.size(); j++) {
            if (tempResult[i].first == result[j].first) {
                result[j].second.numerator =
                    result[j].second.numerator * tempResult[i].second.denominator +
                    tempResult[i].second.numerator * result[j].second.denominator;
                result[j].second.denominator *= tempResult[i].second.denominator;

                int gcd = GCD(abs(result[j].second.numerator), abs(result[j].second.denominator));
                result[j].second.numerator /= gcd;
                result[j].second.denominator /= gcd;

                if (result[j].second.denominator < 0) {
                    result[j].second.numerator = -result[j].second.numerator;
                    result[j].second.denominator = -result[j].second.denominator;
                }

                combined = true;
                break;
            }
        }

        if (!combined) {

            result.emplace_back(tempResult[i]);
        }
    }
    sortVector(result);
}

void EvaluatePolynomial(const vector<pair<int, Rational>>& polynomial, const Rational& input, Rational& result) {
    result = createRational(0, 1);

    for (int i = 0; i < polynomial.size(); ++i) {
        const int power = polynomial[i].first;
        const Rational& coefficient = polynomial[i].second;

        Rational powerResult = createRational(1, 1);
        for (int powerStep = 0; powerStep < power; ++powerStep) {
            powerResult.numerator *= input.numerator;
            powerResult.denominator *= input.denominator;
        }

        Rational termResult = createRational(
            coefficient.numerator * powerResult.numerator,
            coefficient.denominator * powerResult.denominator 
        );

        result.numerator = result.numerator * termResult.denominator + termResult.numerator * result.denominator;
        result.denominator *= termResult.denominator;

        int gcd = GCD(abs(result.numerator), abs(result.denominator));
        result.numerator /= gcd;
        result.denominator /= gcd;

        if (result.denominator < 0) {
            result.numerator = -result.numerator;
            result.denominator = -result.denominator;
        }
    }
}

Rational simplifyFraction(const Rational& r) {
   int gcd = GCD(abs(r.numerator), abs(r.denominator));
   return { r.numerator / gcd, r.denominator / gcd };
}

void DividePolynomials(const vector<pair<int, Rational>>& P1, const vector<pair<int, Rational>>& P2, vector<pair<int, Rational>>& Q, vector<pair<int, Rational>>& R) {
    Q.clear();
    R = P1;

    while (!R.empty() && R.front().first >= P2.front().first) {
        int powerDiff = R.front().first - P2.front().first;
        Rational coeffQuotient = divideRationals(R.front().second, P2.front().second);

        Q.emplace_back(powerDiff, coeffQuotient);

        vector<pair<int, Rational>> tempPoly;
        for (const auto& term : P2) {
            tempPoly.emplace_back(term.first + powerDiff, createRational(term.second.numerator * coeffQuotient.numerator, term.second.denominator * coeffQuotient.denominator));
        }

        vector<pair<int, Rational>> newR;
        size_t i = 0, j = 0;
        while (i < R.size() || j < tempPoly.size()) {
            if (i < R.size() && (j >= tempPoly.size() || R[i].first > tempPoly[j].first)) {
                newR.emplace_back(R[i]);
                i++;
            }
            else if (j < tempPoly.size() && (i >= R.size() || tempPoly[j].first > R[i].first)) {
                newR.emplace_back(tempPoly[j].first, createRational(-tempPoly[j].second.numerator, tempPoly[j].second.denominator));
                j++;
            }
            else {
                Rational diff = {
                    R[i].second.numerator * tempPoly[j].second.denominator - tempPoly[j].second.numerator * R[i].second.denominator,
                    R[i].second.denominator * tempPoly[j].second.denominator };
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

bool getValidatedNumber(const char* input, int& number) {
    for (size_t i = 0; input[i] != '\0'; ++i) {
        if (i >= 2 || input[i] < '0' || input[i] > '9') { 
            cout << "Invalid input."<<endl;
            return false;
        }
    }

    number = 0;
    for (size_t i = 0; input[i] != '\0'; ++i) {
        number = number * 10 + (input[i] - '0');
    }

    return true;
}

void VietaFormulas(const vector<pair<int, Rational>>& polynomial) {
    if (polynomial.empty()) {
        cout << "The polynomial is empty!" << endl;
        return;
    }

    cout << "P(x) = ";
    //DisplayPolynomial(const_cast<vector<pair<int, Rational>>&>(polynomial));

    cout << "Vieta's Formulas for polynomial P(x):" << endl;

    vector<pair<int, Rational>> sortedPolynomial = polynomial;
    sortVector(sortedPolynomial);

    int degree = sortedPolynomial[0].first;    
    Rational leadingCoeff = sortedPolynomial[0].second; 

    for (size_t i = 1; i < sortedPolynomial.size(); ++i) {
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

Rational powRational(const Rational& base, int exponent) {
    if (exponent == 0) return createRational(1, 1);

    Rational result = base;
    for (int i = 1; i < exponent; i++) {
        result = multiplyRationals(result, base);
    }
    return result;
}

Rational BinomialCoefficient(int n, int k) {
    if (k == 0 || k == n) return createRational(1, 1);

    Rational result = createRational(1, 1); 

    for (int i = 1; i <= k; i++) {
        result = multiplyRationals(result, createRational(n - i + 1, 1));
        result = divideRationals(result, createRational(i, 1));
    }

    return result;
}

void RepresentPolynomialInPowersOfXPlusA(vector<pair<int, Rational>>& polynomial, const Rational& shift) {
    if (polynomial.empty()) {
        cout << "The polynomial is empty!" << endl;
        return;
    }

    cout << "P(x) = ";
    DisplayPolynomial(polynomial);
    cout << "Shift value (a): " << shift.numerator << "/" << shift.denominator << endl;

    vector<pair<int, Rational>> expandedPolynomial;

    for (size_t i = 0; i < polynomial.size(); i++) {
        int power = polynomial[i].first;
        Rational coeff = polynomial[i].second;

        for (int k = 0; k <= power; k++) {
            Rational binomialCoeff = BinomialCoefficient(power, k);

            Rational shiftTerm = powRational(shift, power - k);
            if ((power - k) % 2 != 0) {
                shiftTerm.numerator = -shiftTerm.numerator;
            }
            Rational termCoeff = multiplyRationals(binomialCoeff, shiftTerm);
            termCoeff = multiplyRationals(termCoeff, coeff);

            bool termFound = false;
            for (size_t j = 0; j < expandedPolynomial.size(); j++) {
                if (expandedPolynomial[j].first == k) {
                    expandedPolynomial[j].second = addRationals(expandedPolynomial[j].second, termCoeff);
                    termFound = true;
                    break;
                }
            }
            if (!termFound) {
                expandedPolynomial.emplace_back(make_pair(k, termCoeff));
            }
        }
    }

    sortVector(expandedPolynomial);

    cout << "P(x + " << shift.numerator;
    if (shift.denominator != 1) {
        cout << "/" << shift.denominator;
    }
    cout << ") = ";

    bool firstTerm = true;
    for (size_t i = 0; i < polynomial.size(); i++) {
        int power = polynomial[i].first;
        Rational coeff = polynomial[i].second;

        for (int k = 0; k <= power; k++) {
            Rational binomialCoeff = BinomialCoefficient(power, k);

            Rational shiftTerm = powRational(shift, power - k);
            if ((power - k) % 2 != 0) {
                shiftTerm.numerator = -shiftTerm.numerator;
            }
            Rational termCoeff = multiplyRationals(binomialCoeff, shiftTerm);
            termCoeff = multiplyRationals(termCoeff, coeff);

            bool termFound = false;
            for (size_t j = 0; j < expandedPolynomial.size(); j++) {
                if (expandedPolynomial[j].first == k) {
                    expandedPolynomial[j].second = addRationals(expandedPolynomial[j].second, termCoeff);
                    termFound = true;
                    break;
                }
            }
            if (!termFound) {
                expandedPolynomial.emplace_back(make_pair(k, termCoeff));
            }
        }
    }

    cout << endl;
}



int main() {
    cout << "Welcome to the Polynomial Calculator!" << endl;

    int pPowers[MAX_ARR_SIZE] = { 0 };
    int qPowers[MAX_ARR_SIZE] = { 0 };
    char input[MAX_ARR_SIZE];
    int pCounter = 0, qCounter = 0;

    vector<pair<int, Rational>> pVector;
    vector<pair<int, Rational>> qVector;
    vector<pair<int, Rational>> result;
    vector<pair<int, Rational>> quotient;
    vector<pair<int, Rational>> remainder;

    while (true) {
        cout << "\nOptions:" << endl;
        cout << "1 - Sum two polynomials" << endl;
        cout << "2 - Subtract two polynomials" << endl;
        cout << "3 - Multiply polynomial by scalar" << endl;
        cout << "4 - Multiply two polynomials" << endl; 
        cout << "5 - Find  value of polynomial at a given number" << endl;
        cout << "6 - Divide two polynomials" << endl;
        cout << "7. GCD of Polynomials"<<endl;
        cout << "8. Display Vieta's formulas for a given polynomial"<<endl;
        cout << "9. Represent a polynomial in powers of (x+a)"<<endl;
        cout << "11 - Quit the application" << endl;
            
        int choice = 0;
        bool validChoice = false;

        do {
            cout << "Enter your choice: ";
            getInput(input);
            validChoice = getValidatedNumber(input, choice);
        } while (!validChoice);

        switch (choice) {
        case 1:
            cout << "Enter powers for P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Enter powers for Q(x):" << endl;
            EnterPowers(qPowers, input, qCounter);

            EnterCoefficients(qPowers, qCounter, qVector, input);
            sortVector(qVector);

            PolynomSum(&pVector, &qVector, result);
            cout << "Result of P(x) + Q(x): ";
            DisplayPolynomial(result);
            break;

        case 2:
            cout << "Enter powers for P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Enter powers for Q(x):" << endl;
            EnterPowers(qPowers, input, qCounter);

            EnterCoefficients(qPowers, qCounter, qVector, input);
            sortVector(qVector);

            PolynomSubtract(&pVector, &qVector, result);
            cout << "Result of P(x) - Q(x): ";
            DisplayPolynomial(result);
            break;      

        case 3: {
            Rational scalar;
            bool validScalar = false;

            do {
                cout << "Enter scalar (as a rational number): ";
                getInput(input);

                if (isCoefficientValid(input)) {
                    if (ProcessRational(input, scalar)) {
                        validScalar = true;
                    }
                    else {
                        cout << "Invalid scalar. Denominator cannot be zero. Please re-enter." << endl;
                    }
                }
                else {
                    cout << "Invalid scalar format. Please re-enter." << endl;
                }
            } while (!validScalar);

            cout << "Enter powers for the polynomial P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            MultiplyByScalar(scalar, pVector);
            cout << "Result of P(x) multiplied by " << scalar.numerator << "/" << scalar.denominator << ": ";
            DisplayPolynomial(pVector);
            break;
        }

        case 4: { 
            cout << "Enter powers for P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Enter powers for Q(x):" << endl;
            EnterPowers(qPowers, input, qCounter);

            EnterCoefficients(qPowers, qCounter, qVector, input);
            sortVector(qVector);

            MultiplyPolynomials(pVector, qVector, result); 
            cout << "Result of P(x) * Q(x): ";
            DisplayPolynomial(result);
            break;
        }

        case 5: {
            cout << "Enter powers for the polynomial P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            Rational inputValue;
            bool validInput = false;

            do {
                cout << "Enter the value of x (as a rational number): ";
                getInput(input);

                if (isCoefficientValid(input) && ProcessRational(input, inputValue)) {
                    validInput = true;
                }
                else {
                    cout << "Invalid input. Please re-enter." << endl;
                }
            } while (!validInput);

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
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Enter powers for Q(x):" << endl;
            EnterPowers(qPowers, input, qCounter);

            EnterCoefficients(qPowers, qCounter, qVector, input);
            sortVector(qVector);

            vector<pair<int, Rational>> Q, R;
            DividePolynomials(pVector, qVector, Q, R);

            cout << "Quotient: ";
            DisplayPolynomial(Q);

            cout << "Remainder: ";
            DisplayPolynomial(R);
            break;
            
        }
        case 7: {
            cout << "Enter powers for P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Enter powers for Q(x):" << endl;
            EnterPowers(qPowers, input, qCounter);

            EnterCoefficients(qPowers, qCounter, qVector, input);
            sortVector(qVector);

            vector<pair<int, Rational>> gcdResult;
            PolynomialGCD(pVector, qVector, gcdResult);

            cout << "GCD of P(x) and Q(x): ";
            DisplayPolynomial(gcdResult);
            break;
        }
        case 8: {
            cout << "Enter powers for the polynomial P(x):" << endl;

            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            cout << "Calculating Vieta's Formulas for P(x):" << endl;
            VietaFormulas(pVector);
            break;
        }
        case 9: {
            cout << "Enter powers for P(x):" << endl;
            EnterPowers(pPowers, input, pCounter);

            EnterCoefficients(pPowers, pCounter, pVector, input);
            sortVector(pVector);

            Rational shiftValue;
            bool validShift = false;

            do {
                cout << "Enter rational number (shift value a): ";
                getInput(input);

                if (isCoefficientValid(input)) {
                    if (ProcessRational(input, shiftValue)) {
                        validShift = true;
                    }
                    else {
                        cout << "Invalid rational number. Please re-enter." << endl;
                    }
                }
                else {
                    cout << "Invalid format for rational number. Please re-enter." << endl;
                }
            } while (!validShift);

            RepresentPolynomialInPowersOfXPlusA(pVector, shiftValue);
            break;
        }


        case 11:
            cout << "Thank you for using the Polynomial Calculator. Goodbye!" << endl;
            return 0;


        default:
            cout << "Invalid choice. Please enter available ones:" << endl;
            break;
        }

        pVector.clear();
        qVector.clear();
        result.clear();
        quotient.clear();
        remainder.clear();
        pCounter = 0;
        qCounter = 0;
        fill(begin(pPowers), end(pPowers), 0);
        fill(begin(qPowers), end(qPowers), 0);
    }

    return 0;
}