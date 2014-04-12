#pragma once

#include <condition_variable>
#include <vector>
#include <stack>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <thread>

namespace gns {
  const size_t UNDEFINED = std::numeric_limits<size_t>::max();
  const size_t INF = std::numeric_limits<size_t>::max();
  struct UndirectedEdge {
    size_t from, to, weight;
    UndirectedEdge(size_t from, size_t to, size_t weight): from(from), to(to), weight(weight) {}
    UndirectedEdge() {}
  };  
  struct Edge {
    size_t to, weight;
    Edge(size_t to, size_t weight): to(to), weight(weight) {}
    explicit Edge(const UndirectedEdge &edge): to(edge.to), weight(edge.weight) {}
    Edge(): to(0), weight(0) {}
  };
  typedef std::vector < UndirectedEdge > EdgesList;
  typedef std::vector < std::vector < Edge > > AdjencyList;
  class Graph {
    AdjencyList adjacencyList;
  public:
    explicit Graph(size_t vertexNumber): adjacencyList(vertexNumber) {}
    Graph() {}
    void addEdge(const UndirectedEdge &edge) {
      adjacencyList[edge.from].push_back(Edge(edge));
    }
    void addEdge(size_t from, size_t to, size_t weight) {
      adjacencyList[from].push_back(Edge(to, weight));
    }
    std::vector < Edge > &operator[] (size_t vertexNumber){
      return adjacencyList[vertexNumber];
    }
    const std::vector < Edge > &operator[] (size_t vertexNumber) const{
      return adjacencyList[vertexNumber];
    }
    void resize(size_t size) {
      adjacencyList.resize(size);
    }
    size_t size() const {
      return adjacencyList.size();
    }
  };
  void merge(Graph &updatedGraph, const Graph &addedGraph) {
    updatedGraph.resize(std::max(updatedGraph.size(),
                                               addedGraph.size()));
    for(size_t k = 0; k < addedGraph.size(); k++)
      updatedGraph[k].insert(updatedGraph[k].end(), 
        addedGraph[k].begin(), addedGraph[k].end());
  }
  Graph invert(const Graph &graph) {
    Graph invertedGraph(graph.size());
    for(size_t i = 0; i < graph.size(); i++)
      for(size_t j = 0; j < graph[i].size(); j++)
        invertedGraph.addEdge(graph[i][j].to, i, graph[i][j].weight);
    return invertedGraph;
  }

  namespace prns {
    template < class EnterFunc, class VisitFunc, class LeaveFunc >
    void DFS(const Graph &graph, size_t root, EnterFunc onEnter, VisitFunc goByTheEdge, LeaveFunc onLeave) {
      if(onEnter(root))
        return;
      for(size_t i = 0; i < graph[root].size(); i++) {
        size_t to = graph[root][i].to;
        size_t weight = graph[root][i].weight;
        if(goByTheEdge(UndirectedEdge(root, to, weight)))
          DFS(graph, to, onEnter, goByTheEdge, onLeave);
      }
      onLeave(root);
    }
    void componentTreeSearch(const Graph &graph, std::vector < bool > &used,
        Graph &tree, size_t vertex, const std::vector < size_t > &colors, size_t color) {
      DFS(graph, vertex, [&used, &colors, color] (size_t vertex) {
        if(used[vertex] || colors[vertex] != color)
          return true;
        used[vertex] = true;
        return false;
      }, [&used, &colors, &tree, color](UndirectedEdge edge) {
        if(!used[edge.to] && colors[edge.to] == color) {
          tree.addEdge(edge);
          return true;
        }
        return false;
      }, [](size_t vertex) {} );
    }
    Graph subGraphWithZeroEdges(const Graph &graph) {
      Graph ans(graph.size());
      for(size_t i = 0; i < graph.size(); i++)
        for(size_t j = 0; j < graph[i].size(); j++)
          if(graph[i][j].weight == 0)
            ans[i].push_back(graph[i][j]);
      return ans;
    }
    void getDFS_vertexList(const Graph &graph, std::stack < size_t > &out,
        std::vector<bool> &used, size_t vertex = 0) {
      DFS(graph, vertex, [&used](size_t vertex) {
        if(used[vertex])
          return true;
        used[vertex] = true;
        return false;
      }, [](UndirectedEdge edge) { return true; }, 
        [&out](size_t vertex) {
        out.push(vertex);
      });
    }
    void paint(const Graph &graph, std::vector <size_t> &colors, size_t color, size_t vertex) {
      DFS(graph, vertex, [&colors, color](size_t vertex) {
        if(colors[vertex] != UNDEFINED)
          return true;
        colors[vertex] = color;
        return false;
      }, [](UndirectedEdge edje) { return true; }, [](size_t vertex){});
    }
    void restruct(Graph &graph, std::vector < size_t > &edgeMinWeight) {
      size_t size = graph.size();
      edgeMinWeight.assign(size, INF);
      for(size_t i = 0; i < graph.size(); ++i)
        for(size_t j = 0; j < graph[i].size(); j++)
          edgeMinWeight[graph[i][j].to] = std::min(edgeMinWeight[graph[i][j].to], graph[i][j].weight);
      for(size_t i = 0; i < graph.size(); ++i)
        for(size_t j = 0; j < graph[i].size(); j++)
          graph[i][j].weight -= edgeMinWeight[graph[i][j].to];
    }
    void restoreFromRestruct(Graph &graph, const std::vector < size_t > &edgeMinWeight) {
      for(size_t i = 0; i < graph.size(); ++i)
        for(size_t j = 0; j < graph[i].size(); j++)
          graph[i][j].weight += edgeMinWeight[graph[i][j].to];
    }
    void getTree(const Graph &graph, std::vector < bool > &used, Graph &tree, size_t root) {
      DFS(graph, root, [&used](size_t vertex) {
        if(used[vertex])
          return true;
        used[vertex] = true;
        return false;
      }, [&used, &tree](UndirectedEdge edge) {
        if(edge.weight == 0 && !used[edge.to]) {
          tree.addEdge(edge);
          return true;
        }
        return false;
      }, [](size_t vertex) {});
    }
    void buildMSTfromCondence(Graph &cond, Graph &tree, Graph &zero, Graph &ans, 
        std::vector < EdgesList > &mainEdges, std::vector < size_t > colors, size_t color, size_t root) {
      size_t size = zero.size();
      std::vector<bool> isExist(color, false);
      std::vector<bool> used(size, false);
      for(size_t i = 0; i < color; i++) {
        for(size_t j = 0; j < tree[i].size(); j++)
          isExist[tree[i][j].to] = true;
        for(size_t j = 0; j < cond[i].size(); j++) {
          size_t to = cond[i][j].to;
          if(isExist[to]) {
            isExist[to] = false;
            UndirectedEdge edge = mainEdges[i][j];
            ans.addEdge(edge);
            Graph compTree(size);
            componentTreeSearch(zero, used, ans, edge.to, colors, colors[edge.to]);
          }
        }
      }
      Graph compTree(size);
      componentTreeSearch(zero, used, ans, root, colors, colors[root]);
    }
    void threadFord_Bellman(std::pair<size_t, size_t> arg1, 
                            std::pair<AdjencyList *, std::vector<size_t> *> arg2,
                            std::pair<size_t *, std::condition_variable *> arg3,
                            std::pair<size_t, std::mutex *> arg4) {
      size_t begin = arg1.first;
      size_t end = arg1.second;
      AdjencyList *adjacencyList = arg2.first;
      std::vector < size_t > *len = arg2.second;
      size_t *finished = arg3.first; 
      std::condition_variable *condition = arg3.second;
      size_t threadsAmount = arg4.first;
      std::mutex *finishedMutex = arg4.second;
      for(size_t k = 0; k < adjacencyList->size(); k++) {
        for(size_t i = begin; i < end; i++)
          for(size_t j = 0; j < (*adjacencyList)[i].size(); j++)
            if((*len)[(*adjacencyList)[i][j].to] != UNDEFINED)
              (*len)[i] = std::min((*len)[(*adjacencyList)[i][j].to] + (*adjacencyList)[i][j].weight, (*len)[i]);

        std::unique_lock<std::mutex> lock(*finishedMutex); 
        (*finished)++;
        lock.unlock();
        if((*finished) < threadsAmount * (k + 1)) {
          std::mutex my_mutex;
          std::unique_lock<std::mutex> my_lock(my_mutex);
          condition->wait(my_lock, [finished, threadsAmount, k] {
            return (*finished) >= threadsAmount * (k + 1);
          });
        } else {
          condition->notify_all();
        }
      }
    }
  };

  Graph getComponentTree(const Graph &graph, size_t vertex, const std::vector < size_t > &colors) {
    size_t size = graph.size();
    Graph tree(size);
    std::vector < bool > used(size, false);
    prns::componentTreeSearch(graph, used, tree, vertex, colors, colors[vertex]);
    return tree;
  }
  size_t splitIntoComponents(const Graph &graph, std::vector < size_t > &colors) {
    size_t size = graph.size();
    colors.assign(size, UNDEFINED);
    std::vector < bool > used(size, false);
    std::stack < size_t > st;
    for(size_t i = 0; i < graph.size(); i++)
      prns::getDFS_vertexList(graph, st, used, i);
    Graph inverted = invert(graph);
    size_t color = 0;
    while(!st.empty()) {
      size_t next = st.top();
      st.pop();
      if(colors[next] == UNDEFINED) {
        prns::paint(inverted, colors, color++, next);
      }
    }
    return color;
  }
  std::vector < EdgesList > restoreEdgesList(const Graph &graph, Graph &cond,
      const std::vector < size_t > &colors, size_t color) {
    cond = Graph(color);
    std::vector < EdgesList > mainEdges(color);
    Graph fullCond(color);
    std::vector < EdgesList > prevEdges(color);
    for(size_t i = 0; i < graph.size(); i++)
      for(size_t j = 0; j < graph[i].size(); j++) {
       size_t to = graph[i][j].to;
       size_t weight = graph[i][j].weight;
       fullCond.addEdge(colors[i], colors[to], weight);
       prevEdges[colors[i]].push_back(UndirectedEdge(i, to, weight));
     }
    std::vector < Edge > dirLst(color);
    EdgesList undirLst(color);
    std::vector < bool > isExist(color, false);
    for(size_t i = 0; i < color; i++) {
      for(size_t j = 0; j < fullCond[i].size(); j++) {
        size_t to = fullCond[i][j].to;
        size_t weight = fullCond[i][j].weight;
        if(!isExist[to] || dirLst[to].weight > weight) {
        isExist[to] = true;
        dirLst[to] = Edge(to, weight);
        undirLst[to] = prevEdges[i][j];
        }
      }
      for(size_t j = 0; j < fullCond[i].size(); j++) {
        size_t to = fullCond[i][j].to;
        if(isExist[to]) {
          isExist[to] = false;
          if(to != i) {
             cond[i].push_back(dirLst[to]);
            mainEdges[i].push_back(undirLst[to]);
          }
        }
      }
    }
    return mainEdges;
  }
  bool isConnected(const Graph &graph, size_t root) {
    size_t size = graph.size();
    std::vector < size_t > colors(size, UNDEFINED);
    prns::paint(graph, colors, 0, root);
    for(size_t i = 0; i < size; i++)
      if(colors[i])
        return false;
    return true;
  }
  bool isTree(const Graph &graph, size_t root) {
    size_t size = graph.size();
    size_t sum = 0;
    for(size_t i = 0; i < size; i++)
      sum += graph[i].size();
    return (sum == size - 1) && isConnected(graph, root);
  }
  Graph getMST(Graph &graph, size_t root) {
    size_t size = graph.size();
    std::vector < size_t > add(size, 0);
    std::vector < bool > used(size, false);
    prns::restruct(graph, add);
    Graph tree(size), ans(size);
    prns::getTree(graph, used, tree, root);
    if(isTree(tree, root)) {
      prns::restoreFromRestruct(graph, add);
      prns::restoreFromRestruct(tree, add);
      return tree;
    }
    Graph cond;
    std::vector < size_t > colors;
    Graph zero = prns::subGraphWithZeroEdges(graph);
    size_t color = splitIntoComponents(zero, colors);
    std::vector < EdgesList > mainEdges = restoreEdgesList(graph, cond, colors, color);
    tree = getMST(cond, colors[root]);
    prns::buildMSTfromCondence(cond, tree, zero, ans, mainEdges, colors, color, root);
    prns::restoreFromRestruct(graph, add);
    prns::restoreFromRestruct(ans, add);
    return ans;
  }
  std::istream &operator>>(std::istream &in, Graph &graph) {
    size_t n, m;
    in >> n >> m;
    graph = Graph(n);
    for(size_t i = 0; i < m; i++) {
    size_t from, to, weight;
    in >> from >> to >> weight;
    from--; to--;
    graph[from].push_back(Edge(to, weight));
    }
    return in;
  }
  std::ostream &operator<<(std::ostream &out, const Graph &graph) {
    for(size_t i = 0; i < graph.size(); i++) {
    out << i + 1 << ": ";
    for(size_t j = 0; j < graph[i].size(); j++)
      out << graph[i][j].to + 1 << '(' << graph[i][j].weight << ") ";
    out << std::endl;
    }
    return out;
  }
  Graph getRandomGraph(size_t size, size_t minAmount, size_t maxAmount, size_t minWeight, size_t maxWeight) {
    Graph graph(size);
    srand((unsigned int)time(NULL));
    size_t count = minAmount + rand() % (maxAmount - minAmount + 1);
    int deltaW = maxWeight - minWeight + 1;
    for(size_t i = 0; i < count; i++) {
      size_t from = rand() % size;
      size_t to = rand() % size;
      int weight = minWeight + rand() % deltaW;
      graph.addEdge(from, to, weight);
    }
    return graph;
  }
  Graph getRandomTree(size_t minSize, size_t maxSize, size_t minWeight, size_t maxWeight) {
	  srand((unsigned int)time(NULL));
	  size_t size = minSize + rand() % (maxSize - minSize + 1);
    Graph tree(size);
	  int deltaW = maxWeight - minWeight + 1;
	  std::vector < size_t > vertex(size);
    for(size_t i = 0; i < size; i++)
      vertex[i] = i;
    for(size_t firstUnused = 1; firstUnused < size; firstUnused++) {
      size_t nextUsed = rand() % firstUnused;
      size_t nextUnused = firstUnused + rand() % (size - firstUnused);
      size_t weight = minWeight + rand() % deltaW;
      tree.addEdge(vertex[nextUsed], vertex[nextUnused], weight);
      std::swap(vertex[firstUnused], vertex[nextUnused]);
    }
    return tree;
  }
  size_t getTotalWeight(const Graph &graph) {
    size_t weight = 0;
    for(size_t i = 0; i < graph.size(); i++)
      for(size_t j = 0; j < graph[i].size(); j++)
        weight += graph[i][j].weight;
    return weight;
  }
  std::vector < size_t > Ford_Bellman(const Graph &graph, size_t root) {
    size_t size = graph.size();
    std::vector < size_t > len(size, INF);
    len[root] = 0;
    for(size_t i = 0; i < size; i++)
      for(size_t j = 0; j < size; j++)
        if(len[j] != INF)
          for(size_t k = 0; k < graph[j].size(); k++) {
            size_t to = graph[j][k].to;
            size_t weight = graph[j][k].weight;
            len[to] = std::min(len[to], len[j] + weight);
          }
    return len;
  }
  std::vector < size_t > multiThreadedFord_Bellman(const Graph &graph, size_t root) {
    size_t size = graph.size();
    AdjencyList adjacencyList(size);
    for(size_t i = 0; i < size; i++)
      for(size_t j = 0; j < graph[i].size(); j++)
        adjacencyList[graph[i][j].to].push_back(Edge(i, graph[i][j].weight));
    const size_t threadsAmount = 5;
    size_t leftInd[threadsAmount], rightInd[threadsAmount];
    for(size_t i = 0; i < threadsAmount; i++)
      leftInd[i] = i * size / threadsAmount;
    for(size_t i = 0; i < threadsAmount - 1; i++)
      rightInd[i] = leftInd[i + 1];
    rightInd[threadsAmount - 1] = size;
    std::vector < std::thread* > threads(threadsAmount);
    std::vector < size_t > len(size, UNDEFINED);
    len[root] = 0;
    std::mutex finishedMutex;
    std::condition_variable condition;
    size_t finished = 0;
    for(size_t j = 0; j < threadsAmount; j++)
      threads[j] = new std::thread(prns::threadFord_Bellman, 
                                   std::make_pair(leftInd[j], rightInd[j]),
                                   std::make_pair(&adjacencyList, &len), 
                                   std::make_pair(&finished, &condition),
                                   std::make_pair(threadsAmount, &finishedMutex));
    for(size_t j = 0; j < threadsAmount; j++) {
      threads[j]->join();
      delete threads[j];
    }
    return len;
  }
};