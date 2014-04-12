#include <iostream>
#include <fstream>
#include "test.h"

using namespace std;

int main() {
  gns::Graph g;
  ifstream in("test_MST.txt");
  ofstream out("output.txt");
  //testGetMST(200);
  testMultiThreadedFord_Bellman(7);
  system("pause");
  in >> g;
  size_t root;
  in >> root;
  vector < size_t > len = multiThreadedFord_Bellman(g, root - 1);
  for(size_t i = 0; i < len.size(); i++)
    out << len[i] << ' ';
  //graph tree;
  //getBrutMST(g, root - 1, tree);
  //out << tree;
  //out << getMST(g, root - 1);
  return 0;
}