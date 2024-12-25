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

void sortVector(vector<pair<Rational, int>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec.size() - i - 1; ++j) {
            if (vec[j].second < vec[j + 1].second) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

bool isInputValid(const char input[]) {
    if (input[0] == '\0') {
        cout << "Input is empty!" << endl;
        return false;
    }

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
    if (!isInputValid(input)) {
        return false;
    }

    for (int i = 0; input[i] != '\0'; ++i) {
        char c = input[i];

        if (c == ' ') continue;

        if (c == '-' || c == '/') {
            cout << "Input contains invalid power format (negative or fractional)!" << endl;
            return false;
        }
    }

    return true;
}

bool isCoefficientValid(const char input[]) {
    if (!isInputValid(input)) {
        return false;
    }

    bool hasSlash = false;
    bool hasMinus = false;
    bool inNumber = false;
    int length = 0;

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
            inNumber = true;
        }
        else if (c == ' ') {
            if (inNumber) {
                cout << "Invalid coefficient: Spaces are not allowed within a number!" << endl;
                return false;
            }
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

    const char* slashPos = strchr(input, '/');
    if (slashPos != nullptr && atoi(slashPos + 1) == 0) {
        cout << "Invalid coefficient: Denominator cannot be zero!" << endl;
        return false;
    }

    return true;
}

int splitInput(const char* input, int* output) {
    if (input == nullptr || output == nullptr) {
        return 0;
    }

    int index = 0;       
    bool inNum = false;  
    bool isNeg = false; 
    int num = 0;         
    bool hasMinus = false;

    while (*input != '\0') {
        if (*input == '-') {
            if (inNum) {
                cout << "Invalid input: multiple '-' in a single number" <<endl;
                return -1;
            }
            isNeg = true;
            hasMinus = true;
            inNum = true;
        }
        else if (*input >= '0' && *input <= '9') {
            num = num * 10 + (*input - '0');
            inNum = true;
        }
        else if (*input == ' ') {
            if (inNum) {
                if (isNeg) {
                    num = -num;
                }
                output[index++] = num; 
                num = 0;               
                inNum = false;
                isNeg = false;
                hasMinus = false;
            }
        }
        input++;
    }

    if (inNum) {
        if (isNeg) {
            num = -num;
        }
        output[index++] = num;
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
            cin.getline(input, 50);

            if (isCoefficientValid(input)) {
                Rational coefficient;
                if (ProcessRational(input, coefficient)) {
                    coefficients.push_back({ powers[i], coefficient });
                    valid = true;
                }
                else {
                    cout << "Invalid coefficient: Denominator cannot be zero. Please re-enter." << endl;
                }
            }
            else {
                cout << "Invalid coefficient. Please re-enter." << endl;
            }
        } while (!valid);
    }
}


void EnterPowers(int* powers, char* input, int& counter) {
    if (powers == nullptr || input == nullptr) {
        cerr << "Error: Null pointer passed to EnterPowers." << endl;
        return;
    }

    do {
        cout << "Enter powers: ";
        cin.getline(input, MAX_ARR_SIZE);
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

int main()
{
    int pPowers[MAX_ARR_SIZE] = { 0 };
    int qPowers[MAX_ARR_SIZE] = { 0 };
    char input[MAX_ARR_SIZE];
    int pCounter = 0;
    int qCounter = 0;

    vector < pair<int, Rational>> pVector;
    vector<pair<int, Rational>> qVector;

    cout << "Enter powers for P(x):" << endl;
    EnterPowers(pPowers, input, pCounter);

    EnterCoefficients(pPowers, pCounter, pVector, input);

    cout << "Enter powers for Q(x):" << endl;
    EnterPowers(qPowers, input, qCounter);

    EnterCoefficients(qPowers, qCounter, qVector, input);


    cout << "P(x):" << endl;
    for (const auto& term : pVector) {
        cout << "Power: " << term.first
            << ", Coefficient: " << term.second.numerator
            << "/" << term.second.denominator << endl;
    }

    
    cout << "Q(x):" << endl;
    for (const auto& term : qVector) {
        cout << "Power: " << term.first
            << ", Coefficient: " << term.second.numerator
            << "/" << term.second.denominator << endl;
    }

    return 0;
}
