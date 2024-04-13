#include <iostream>
#include <algorithm>
#include <string>
#include <limits>
#include <vector>
#include "needleman-wunsch.hpp"

using namespace std;

vector<pair<int, pair<string, string>>> StarAlignment(vector<string> sequences) {
  int m = sequences.size();

  vector<vector<pair<int, pair<string, string>>>> alignments(m, vector<pair<int, pair<string, string>>>(m, make_pair(0, make_pair(" ", " "))));


  for (int i = 0; i < m; i++) {
    for (int j = i + 1; j < m; j++) { 
      pair<int, pair<string, string>> alignment = NeedlemanWunsch(sequences[i], sequences[j]);
      alignments[i][j] = alignment;
      alignments[j][i] = alignment;
    }
  }

  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      cout << alignments[i][j].first << "\t";
    }
    cout << endl;
  }
  
  vector<int> total_scores(m, 0);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      if (i != j) {
        total_scores[i] += alignments[i][j].first;
      }
    }
  }

  int max_index = max_element(total_scores.begin(), total_scores.end()) - total_scores.begin();

  vector<pair<int, pair<string, string>>> final_alignments(m, {0, {" ", " "}});

  for (int i = 0; i < m; i++) {
    if (i != max_index) {
      final_alignments[i] = alignments[i][max_index];
    } else {
      final_alignments[i].first = 0; 
    }
  }


  int max_length = final_alignments[max_index].second.first.length();


  for (int i = 0; i < m; i++) {
      string aligned_seq = final_alignments[i].second.first;
      string gaps(max_length - aligned_seq.length(), '-');
      cout << aligned_seq << gaps << endl;
  }
  
  return final_alignments;
}

int main() {
  vector<string> seqs = read_sequences("B3.txt");
  for (auto& seq : seqs)
    cout << seq << endl;

  vector<pair<int, pair<string, string>>> alignments = StarAlignment(seqs);

}


// 