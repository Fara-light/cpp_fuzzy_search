#include <iostream>
#include <initializer_list>
#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

long long getGapPenalty(long long gapLength) {
  return 2 * gapLength;
}

long long getSimilarityScore(char lhs, char rhs) {
  if (lhs == rhs) {
    return 3;
  }
  return -3;
}

long long getMaxVerticalScore(const std::vector<std::vector<long long>>& substitutionMatrix, size_t i, size_t j) {
  return substitutionMatrix[i - 1][j] - getGapPenalty(1);
}

long long getMaxHorizontalScore(const std::vector<std::vector<long long>>& substitutionMatrix, size_t i, size_t j) {
  return substitutionMatrix[i][j - 1] - getGapPenalty(1);
}

std::vector<std::vector<long long>> getSubstitutionMatrix(std::string_view sequence1, std::string_view sequence2) {
  std::vector<std::vector<long long>> substitutionMatrix(sequence1.length() + 1, std::vector<long long>(sequence2.length() + 1, 0));
  for (size_t i = 1; i < substitutionMatrix.size(); ++i) {
    for (size_t j = 1; j < substitutionMatrix[i].size(); ++j) {
      substitutionMatrix[i][j] = std::max({
        substitutionMatrix[i - 1][j - 1] + getSimilarityScore(sequence1[i - 1], sequence2[j - 1]),
        getMaxHorizontalScore(substitutionMatrix, i, j),
        getMaxVerticalScore(substitutionMatrix, i, j),
        0ll
      });
    }
  }
  return substitutionMatrix;
}

long long getMatrixMaximumValue(const std::vector<std::vector<long long>>& matrix, long long minimumValue) {
  long long maximumValue = minimumValue;
  for (const auto& line: matrix) {
    for (const auto& element: line) {
      if (element > maximumValue) {
        maximumValue = element;
      }
    }
  }
  return maximumValue;
}

std::vector<std::pair<size_t, size_t>> getMatrixMaximumValueCoordinates(const std::vector<std::vector<long long>>& matrix, long long minimumValue) {
  long long maximiumScore = getMatrixMaximumValue(matrix, minimumValue);
  std::vector<std::pair<size_t, size_t>> maximumScoreCoordinates;
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j) {
      if (matrix[i][j] == maximiumScore) {
        maximumScoreCoordinates.emplace_back(i, j);
      }
    }
  }
  return maximumScoreCoordinates;
}

struct SubsequenceInformation {
  SubsequenceInformation(size_t pos, size_t len) : pos(pos), len(len) {}
  size_t pos;
  size_t len;
};

SubsequenceInformation traceFromCoordinates(const std::vector<std::vector<long long>>& substitutionMatrix, size_t i, size_t j) {
  size_t end = i;
  while (substitutionMatrix[i][j] != 0) {
    long long diagonal = substitutionMatrix[i - 1][j - 1];
    long long vertical = substitutionMatrix[i - 1][j];
    long long horizontal = substitutionMatrix[i][j - 1];
    if (diagonal >= vertical && diagonal >= horizontal) {
      --i;
      --j;
    } else if (vertical >= horizontal) {
      --i;
    } else {
      --j;
    }
  }
  if (end == i) {
    return SubsequenceInformation(0, 0);
  }
  return SubsequenceInformation(i, end - i);
}

std::vector<SubsequenceInformation> getSubsequenceInformation(const std::vector<std::vector<long long>>& substitutionMatrix) {
  auto maximumScoreCoordinates = getMatrixMaximumValueCoordinates(substitutionMatrix, 0);
  std::vector<SubsequenceInformation> subsequenceInformation;
  for (const auto& [i, j]: maximumScoreCoordinates) {
    subsequenceInformation.push_back(traceFromCoordinates(substitutionMatrix, i, j));
  }
  return subsequenceInformation;
}

long long calculateEditDistance(std::string_view lhs, std::string_view rhs) {
  std::vector<std::vector<long long>> distanceMatrix(lhs.size() + 1, std::vector<long long>(rhs.size() + 1));
  for (size_t i = 1; i < distanceMatrix.size(); ++i) {
    distanceMatrix[i][0] = i;
  }
  for (size_t j = 1; j < distanceMatrix[0].size(); ++j) {
    distanceMatrix[0][j] = j;
  }
  for (size_t i = 1; i < distanceMatrix.size(); ++i) {
    for (size_t j = 1; j < distanceMatrix[i].size(); ++j) {
      long long substitutionCost = 1;
      if (lhs[i - 1] == rhs[j - 1]) {
        substitutionCost = 0;
      }
      distanceMatrix[i][j] = std::min({
        distanceMatrix[i - 1][j] + 1,
        distanceMatrix[i][j - 1] + 1,
        distanceMatrix[i - 1][j - 1] + substitutionCost
      });
    }
  }
  return distanceMatrix.back().back();
}

long long calculateTargetBestEditDistance(std::string_view lhs, std::string_view rhs) {
  auto substitutionMatrix = getSubstitutionMatrix(lhs, rhs);
  auto subsequenceInformationVector = getSubsequenceInformation(substitutionMatrix);
  long long minimumSubsequenceDistance = calculateEditDistance(lhs.substr(subsequenceInformationVector[0].pos, subsequenceInformationVector[0].len), rhs);
  for (size_t i = 1; i < subsequenceInformationVector.size(); ++i) {
    SubsequenceInformation& subsequenceInformation = subsequenceInformationVector[i];
    std::string_view subsequence = lhs.substr(subsequenceInformation.pos, subsequenceInformation.len);
    long long distance = calculateEditDistance(subsequence, rhs);
    if (distance < minimumSubsequenceDistance) {
      minimumSubsequenceDistance = distance;
    }
  }
  return minimumSubsequenceDistance;
}

std::vector<std::pair<std::string, size_t>> fuzzySortStrings(const std::vector<std::string>& targets, const std::string& search) {
  std::vector<std::pair<std::string, size_t>> targetsEvaluation;
  std::transform(
    targets.begin(),
    targets.end(),
    std::back_inserter(targetsEvaluation),
    [search](const std::string& str) {
      return std::make_pair(str, calculateTargetBestEditDistance(str, search));
    }
  );
  std::sort(
    targetsEvaluation.begin(),
    targetsEvaluation.end(),
    [](auto lhs, auto rhs) {
      return lhs.second < rhs.second;
    }
  );
  return targetsEvaluation;
}

int main() {
  size_t stringNumber;
  std::cout << "Enter string number: ";
  std::cin >> stringNumber;
  std::vector<std::string> targetStrings(stringNumber);
  for (auto& str: targetStrings) {
    std::cin >> str;
  }
  std::string search;
  std::cout << "Enter search query: ";
  std::cin >> search;
  auto searchResult = fuzzySortStrings(targetStrings, search);
  for (const auto& [str, score]: searchResult) {
    std::cout << str << '\t' << score << '\n';
  }
  return 0;
}
