#include <iostream>
#include <vector>

using namespace std;

const int MAX_ARR_SIZE = 50;

struct Rational {
    int numerator;
    int denominator;

    Rational(int num = 0, int den = 1) {
        numerator = num;
        denominator = den;
    }
};

void sortVector(vector<pair<int, Rational>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec.size() - i - 1; ++j) {
            if (vec[j].first < vec[j + 1].first) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

void getInput(char* input) {
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


void DisplayPolynom(vector<pair<int, Rational>>* vPtr) {
    if (vPtr == nullptr) return;
    vector<pair<int, Rational>>& vec = *vPtr; 

    for (size_t i = 0; i < vec.size(); ++i) {
        pair<int, Rational> term = vec[i];  
        int power = term.first;
        Rational coeff = term.second;       

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

bool isInputValid(const char input[]) {
    if ((input == nullptr) ? (cout << "Error: Null pointer passed as input!" << endl, true) : false) return false;

    if ((input[0] == '\0') ? (cout << "Input is empty!" << endl, true) : false) return false;

    bool hasNonSpaceCharacter = false;

    for (int i = 0; input[i] != '\0'; ++i) {
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

    for (int i = 0; input[i] != '\0'; ++i) {
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
    if (!isInputValid(input)) {
        return false;
    }
    bool hasSlash = false, hasMinus = false, inNumber = false;
    int length = 0, denominator = 0;

    for (int i = 0; input[i] != '\0'; ++i, ++length) {
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

    for (int i = 0; input[i] != '\0'; ++i) {
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
            }
            else {
                numerator = numerator * 10 + (input[i] - '0');
            }
        }
        else {
            return false;
        }
    }

    if (inDenominator && denominator == 0) {
        return false; 
    }

    if (isNegative) {
        numerator = -numerator;
    }

    rational = Rational(numerator, inDenominator ? denominator : 1);
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

    return Rational(sumNumerator / gcd, lcm / gcd);
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

    vector<pair<int, Rational>>& vec1 = *v1Ptr, & vec2 = *v2Ptr;

    size_t idV1 = 0, idV2 = 0;

    while (idV1 < vec1.size() && idV2 < vec2.size()) {
        pair<int, Rational> term1 = vec1[idV1], term2 = vec2[idV2];

        if (term1.first == term2.first) {
            Rational diffCoeff = AddRational(term1.second, Rational(-term2.second.numerator, term2.second.denominator));
            resVec.emplace_back(term1.first, diffCoeff);
            idV1++;
            idV2++;
        }
        else if (term1.first > term2.first) { 
            resVec.emplace_back(term1.first, term1.second);
            idV1++;
        }
        else { 
            Rational negCoeff(-term2.second.numerator, term2.second.denominator);
            resVec.emplace_back(term2.first, negCoeff);
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
        Rational negCoeff(-remainingTerm.second.numerator, remainingTerm.second.denominator);
        resVec.emplace_back(remainingTerm.first, negCoeff);
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
void MultiplyPolynomials(
    const vector<pair<int, Rational>>& poly1,
    const vector<pair<int, Rational>>& poly2,
    vector<pair<int, Rational>>& result)
{
    result.clear();

    if (poly1.empty() || poly2.empty()) {
        return;
    }

    vector<pair<int, Rational>> tempResult;

    for (size_t i = 0; i < poly1.size(); ++i) {
        for (size_t j = 0; j < poly2.size(); ++j) {
            Rational coefficient(
                poly1[i].second.numerator * poly2[j].second.numerator,
                poly1[i].second.denominator * poly2[j].second.denominator
            );

            int gcd = GCD(abs(coefficient.numerator), abs(coefficient.denominator));
            coefficient.numerator /= gcd;
            coefficient.denominator /= gcd;

            if (coefficient.denominator < 0) {
                coefficient.numerator = -coefficient.numerator;
                coefficient.denominator = -coefficient.denominator;
            }

            int power = poly1[i].first + poly2[j].first;

            tempResult.push_back({ power, coefficient });
        }
    }

    for (size_t i = 0; i < tempResult.size(); ++i) {
        bool combined = false;
        for (size_t j = 0; j < result.size(); ++j) {
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

            result.push_back(tempResult[i]);
        }
    }
    sortVector(result);
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


int main() {
    cout << "Welcome to the Polynomial Calculator!" << endl;

    int pPowers[MAX_ARR_SIZE] = { 0 };
    int qPowers[MAX_ARR_SIZE] = { 0 };
    char input[MAX_ARR_SIZE];
    int pCounter = 0, qCounter = 0;

    vector<pair<int, Rational>> pVector;
    vector<pair<int, Rational>> qVector;
    vector<pair<int, Rational>> result;

    while (true) {
        cout << "\nOptions:" << endl;
        cout << "1 - Sum two polynomials" << endl;
        cout << "2 - Subtract two polynomials" << endl;
        cout << "3 - Quit the application" << endl;
        cout << "4 - Multiply polynomial by scalar" << endl;
        cout << "5 - Multiply two polynomials" << endl; 

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
            DisplayPolynom(&result);
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
            DisplayPolynom(&result);
            break;

        case 3:
            cout << "Thank you for using the Polynomial Calculator. Goodbye!" << endl;
            return 0;

        case 4: {
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
            DisplayPolynom(&pVector);
            break;
        }

        case 5: { 
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
            DisplayPolynom(&result);
            break;
        }

        default:
            cout << "Invalid choice. Please enter 1, 2, 3, 4, or 5." << endl;
            break;
        }

        pVector.clear();
        qVector.clear();
        result.clear();
        pCounter = 0;
        qCounter = 0;
        fill(begin(pPowers), end(pPowers), 0);
        fill(begin(qPowers), end(qPowers), 0);
    }



    return 0;
}



