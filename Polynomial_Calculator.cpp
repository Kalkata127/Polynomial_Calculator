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

float processFractionInput(char* ptr) {
    if (ptr == nullptr) {
        cout << "ptr is null";
        return 0;
    }
    char* originalPtr = ptr;
    float result = 0.0;
    int slashid = 0;
    int counter = 0;
    bool IsNegatie = 0;

    if (*ptr == '-') {
        ptr++;
        originalPtr++;
        IsNegatie = 1;
    }

    while (*ptr != '\0') {
        if (*ptr == '/')
            slashid = counter;
        ptr++;
        counter++;
    }
    ptr = originalPtr;

    float first = 0.0;
    for (int i = slashid; i > 0; i--) {
        if (i > 1) {
            first += (((float)*ptr - '0') * 10);
        }
        else {
            first += ((float)*ptr - '0');
        }
        ptr++;
    }

    ptr = originalPtr + slashid + 1;

    float second = 0.0;
    for (int j = counter - slashid - 1; j > 0; j--) {
        if (j > 1) {
            second += (((float)*ptr - '0') * 10);
        }
        else {
            second += ((float)*ptr - '0');
        }
        ptr++;
    }
    if (second == 0.00) {
        cout << "Can't divide by 0";
        return 0;
    }
    result = first / second;

    if (IsNegatie) {
        return result * (-1);
    }
    return result;
}

void toInt(char* input, int* output, int& counter) {
    if (input == nullptr || output == nullptr) {
        return;
    }

    bool isNeg = 0;
    int num = 0;
    bool inNumber = 0;

    while (*input != '\0') {
        if (*input == '-') {
            isNeg = true;
        }
        else if (*input >= '0' && *input <= '9') {
            num = num * 10 + (*input - '0');
            inNumber = true;
        }
        else if (*input == ' ') {
            if (inNumber) {
                if (isNeg) {
                    num = -num;
                }
                *output = num;
                output++;
                counter++;
                num = 0;
                isNeg = false;
                inNumber = false;
            }
        }
        input++;
    }

    if (inNumber) {
        if (isNeg) {
            num = -num;
        }
        *output = num;
        counter++;
    }
}

void DisplayPolynom(vector<pair<int, int>>* vPtr) {
    if (vPtr == nullptr) return;

    vector<pair<int, int>>* orgPtr = vPtr;
    for (size_t i = 0; i < (*vPtr).size(); i++) {
        if ((*vPtr)[i].second > 0) {
            cout << "+" << (*vPtr)[i].second;
        }
        else {
            cout << (*vPtr)[i].second;
        }
        cout << "x^" << (*vPtr)[i].first << " ";
    }
    cout << endl;

}

void sortVector(vector<pair<int, int>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec.size() - i - 1; ++j) {
            if (vec[j].first < vec[j + 1].first) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

void PolynomSum(vector<pair<int, int>>* v1Ptr, vector<pair<int, int>>* v2Ptr, vector<pair<int, int>>& resVec) {
    if (v1Ptr == nullptr || v2Ptr == nullptr) {
        return;
    }
    size_t idV1 = 0;
    size_t idV2 = 0;

    while (idV1 < (*v1Ptr).size() && idV2 < (*v2Ptr).size()) {
        if ((*v1Ptr)[idV1].first == (*v2Ptr)[idV2].first) {
            int sum = (*v1Ptr)[idV1].second + (*v2Ptr)[idV2].second;
            resVec.emplace_back((*v1Ptr)[idV1].first, sum);
            idV1++;
            idV2++;
        }
        else if ((*v1Ptr)[idV1].first > (*v2Ptr)[idV2].first) {
            resVec.emplace_back((*v1Ptr)[idV1].first, (*v1Ptr)[idV1].second);
            idV1++;
        }
        else {
            resVec.emplace_back((*v2Ptr)[idV2].first, (*v2Ptr)[idV2].second);
            idV2++;
        }
    }

    while (idV1 < (*v1Ptr).size()) {
        resVec.emplace_back((*v1Ptr)[idV1].first, (*v1Ptr)[idV1].second);
        idV1++;
    }

    while (idV2 < (*v2Ptr).size()) {
        resVec.emplace_back((*v2Ptr)[idV2].first, (*v2Ptr)[idV2].second);
        idV2++;
    }
}

void PolynomSubtract(vector<pair<int, int>>* v1Ptr, vector<pair<int, int>>* v2Ptr, vector<pair<int, int>>& resVec) {
    if (v1Ptr == nullptr || v2Ptr == nullptr) {
        return;
    }
    size_t idV1 = 0;
    size_t idV2 = 0;

    while (idV1 < (*v1Ptr).size() && idV2 < (*v2Ptr).size()) {
        if ((*v1Ptr)[idV1].first == (*v2Ptr)[idV2].first) {
            int sum = (*v1Ptr)[idV1].second - (*v2Ptr)[idV2].second;
            resVec.emplace_back((*v1Ptr)[idV1].first, sum);
            idV1++;
            idV2++;
        }
        else if ((*v1Ptr)[idV1].first > (*v2Ptr)[idV2].first) {
            resVec.emplace_back((*v1Ptr)[idV1].first, (*v1Ptr)[idV1].second);
            idV1++;
        }
        else {
            resVec.emplace_back((*v2Ptr)[idV2].first, (*v2Ptr)[idV2].second);
            idV2++;
        }
    }

    while (idV1 < (*v1Ptr).size()) {
        resVec.emplace_back((*v1Ptr)[idV1].first, (*v1Ptr)[idV1].second);
        idV1++;
    }

    while (idV2 < (*v2Ptr).size()) {
        resVec.emplace_back((*v2Ptr)[idV2].first, (*v2Ptr)[idV2].second);
        idV2++;
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

int main()
{
    int pPowers[MAX_ARR_SIZE] = { 0 };
    int qPowers[MAX_ARR_SIZE] = { 0 };
    char input[MAX_ARR_SIZE];
    int pCounter = 0;

    do {
        cout << "Enter powers of P(x): ";
        cin.getline(input, MAX_ARR_SIZE);
        pCounter = 0;
        //toInt(input, pPowers, pCounter);
    } while (!arePowersValid(input));

    for (int i = 0; i < pCounter; i++) {
        cout << pPowers[i];
    }

    cout << "Enter powers of Q(x): ";
    cin.getline(input, MAX_ARR_SIZE);
    int qCounter = 0;
    toInt(input, qPowers, qCounter);

    vector<pair<int, int>> v1;

    int i = 0;
    do {
        cout << "Enter P coeficient before " << "x^" << pPowers[i] << " :";
        int psecond = 0;
        cin >> psecond;
        v1.emplace_back(pPowers[i], psecond);
        i++;
    } while (i < pCounter);

    sortVector(v1);

    for (size_t i = 0; i < v1.size(); i++)
        cout << v1[i].first << " " << v1[i].second << endl;

    vector<pair<int, int>> v2;
    int j = 0;
    do {
        cout << "Enter Q coeficient before " << "x^" << qPowers[j] << " :";
        int qsecond = 0;
        cin >> qsecond;
        v2.emplace_back(qPowers[j], qsecond);
        j++;
    } while (j < qCounter);

    sortVector(v2);

    for (size_t i = 0; i < v2.size(); i++)
        cout << v2[i].first << " " << v2[i].second << endl;
    cout << "P(x) = ";
    DisplayPolynom(&v1);
    cout << "Q(x) = ";
    DisplayPolynom(&v2);

    cout << "Sum is: ";
    vector<pair<int, int>> resVec;
    PolynomSum(&v1, &v2, resVec);

    vector<pair<int, int>> resVecSubt;
    PolynomSubtract(&v1, &v2, resVecSubt);
    cout << "Sum: " << endl;
    DisplayPolynom(&resVec);
    cout << endl << "Subtract: " << endl;
    DisplayPolynom(&resVecSubt);
    return 0;
}
