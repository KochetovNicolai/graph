#ifndef Graph_TEST
#define Graph_TEST

#include <ctime>
#include "Graph.h"
using namespace gns;

size_t getBrutMST(const Graph &graph, size_t root, Graph &ans) {
  size_t adjacencyListAmount = 0;
  for(size_t i = 0; i < graph.size(); i++)
    adjacencyListAmount += graph[i].size();
  size_t mask = ((size_t) 1) << adjacencyListAmount;
  size_t size = graph.size();
  for(size_t i = 0; i < mask; i++) {
    size_t sum = 0;
    for(size_t j = 1; j < mask; j <<= 1)
      sum += i & j ? 1 : 0;
    if(sum != size - 1)
      continue;
    Graph tree(size);
    size_t pts = 1;
    for(size_t j = 0; j < size; j++)
      for(size_t k = 0; k < graph[j].size(); k++) {
        if(pts & i)
          tree[j].push_back(graph[j][k]);
        pts <<= 1;
      }
    if(isTree(tree, root) && (ans.size() == 0 || getTotalWeight(tree) < getTotalWeight(ans)))
      ans = tree;
  } 
  return getTotalWeight(ans);
}
bool testGetMST(size_t testAmount) {
  char pref[] = "MST test: ";
  for(size_t i = 0; i < testAmount; i++) {
    Graph tree = getRandomTree(8, 12, 1, 20);
    if(!isTree(tree, 0)) {
      std::cout << pref << "Getting random tree failed on test " << i + 1 << std::endl;
      return false;
    }
    merge(tree, getRandomGraph(tree.size(), 0, 10, 1, 20));
    Graph tree1, tree2;
    tree1 = getMST(tree, 0);
    if(getBrutMST(tree, 0, tree2) != getTotalWeight(tree1)) {
      std::cout << pref << "Random test " << i + 1 << " faled!" << std::endl;
      std::cout << pref << "Brutforse MST:\n" << tree2 << "\nAlgo MST:\n" << tree1;
      return false;
    }
    std::cout << pref << "test " << i + 1 << " OK" << std::endl;
  }
  return true;
}
template < class T >
std::ostream &operator << (std::ostream &out, const std::vector < T > &v) {
  for(size_t i = 0; i < v.size(); i++)
    out << v[i] << ' ';
  return out;
}
bool testMultiThreadedFord_Bellman(size_t testAmount) {
  char pref[] = "Multi-threaded Ford-Bellman test: ";
  double summaryTime = 0;
  for(size_t i = 0; i < testAmount; i++) {
    Graph tree = getRandomTree(1000, 1000, 1, 100);
    if(!isTree(tree, 0)) {
      std::cout << pref << "Getting random tree failed on test " << i + 1 << std::endl;
      return false;
    }
    merge(tree, getRandomGraph(tree.size(), 1000000, 1000000, 1, 100));

    double baseTime = clock();
    std::vector < size_t > len1 = Ford_Bellman(tree, 0);
    double brutTime = clock() - baseTime;
    std::vector < size_t > len2 = multiThreadedFord_Bellman(tree, 0);
    double multiThreadedTime = clock() - brutTime - baseTime;
    summaryTime += brutTime - multiThreadedTime;
    if(len1 != len2) {
      std::cout << pref << "Random test " << i + 1 << " faled!" << std::endl;
      std::cout << pref << "Standart Ford-Bellman:\n" << len1;
      std::cout << pref << "\nMulti-threaded F-B:\n" << len2;
      std::cout << std::endl;
      return false;
    }
    brutTime /= CLOCKS_PER_SEC;
    multiThreadedTime /= CLOCKS_PER_SEC;
    std::cout << pref <<  "test " << i + 1 << " OK " ;
    std::cout << multiThreadedTime << " vs " << brutTime;
    std::cout << " (" << brutTime - multiThreadedTime << " seconds)" << std::endl;
  }
  summaryTime /= CLOCKS_PER_SEC * testAmount;
  std::cout << "---- Average difference in " << summaryTime << " seconds" << std::endl;
  return true;
}
#endif //Graph_TEST