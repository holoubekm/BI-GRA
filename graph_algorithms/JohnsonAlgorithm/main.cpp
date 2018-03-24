#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <memory>
#include <algorithm>
#include <queue>
#include <vector>
#include <iomanip>
#include <stack>
#include <ctime>
#include <set>

using namespace std;

template<typename V>
class Value
{
public:
	Value() : val(NONE)
	{
	}

	Value(V inVal) : val(inVal)
	{
	}

	bool operator ==(const Value<V>& rhs) const
	{
		return val == rhs.val;
	}

	bool operator !=(const Value<V>& rhs) const
	{
		return !operator==(rhs);
	}

	bool operator >(const Value<V>& rhs) const
	{
		return val > rhs.val;
	}

	template<V>
	friend Value<V> operator +(const Value<V>& lhs, const Value<V>& rhs);
	template<V>
	friend Value<V> operator -(const Value<V>& lhs, const Value<V>& rhs);

public:
	V val;

	static const int32_t MAX_VAL = 65534;
	static const int32_t MIN_VAL = -65534;
	static const int32_t NONE = -100000;
};

template<typename Value>
struct Value_Cmp
{
public:
	bool operator()(const Value& a, const Value& b) const
	{
		return a.val < b.val;
	}
};

template<typename V>
Value<V> operator +(const Value<V>& lhs, const Value<V>& rhs)
{
	if (lhs.val == Value<V>::MAX_VAL || rhs.val == Value<V>::MAX_VAL)
		return Value<V>::MAX_VAL;
	if (lhs.val == Value<V>::MIN_VAL || rhs.val == Value<V>::MIN_VAL)
		return Value<V>::MIN_VAL;
	return lhs.val + rhs.val;
}

template<typename V>
Value<V> operator -(const Value<V>& lhs, const Value<V>& rhs)
{
	if (lhs.val == Value<V>::MAX_VAL || rhs.val == Value<V>::MIN_VAL)
		return Value<V>::MAX_VAL;
	if (lhs.val == Value<V>::MIN_VAL || rhs.val == Value<V>::MAX_VAL)
		return Value<V>::MIN_VAL;
	return lhs.val - rhs.val;
}

template<typename T>
class Node
{
public:
	Node(int32_t inx) : state(State::FRESH), dist(T::MAX_VAL), index(inx), parent(NO_PARENT), oTime(NO_TIME), cTime(NO_TIME)
	{
	}

	void setDistance(T inDist)
	{
		dist = inDist;
	}

	void setParent(Node<T>& father)
	{
		parent = father.index;
	}

	void setParent(int32_t father)
	{
		parent = father;
	}

	void close()
	{
		state = State::CLOSED;
	}

	void setCloseTime(int32_t time)
	{
		cTime = time;
	}

	void setOpenTime(int32_t time)
	{
		oTime = time;
	}

	void reset()
	{
		state = State::FRESH;
		parent = NO_PARENT;
		dist = T::MAX_VAL;
		cTime = oTime = NO_TIME;
	}

	operator int()
	{
		return index;
	}

	enum class State
	{
		FRESH,
		OPEN,
		CLOSED
	};
	
public:
	using Value_Type = typename T;
	static const int32_t NO_PARENT = -1;
	static const int32_t NO_TIME = -1;

	State state;
	T dist;
	int32_t index;
	int32_t parent;
	int32_t oTime;
	int32_t cTime;
};

template<typename N>
struct Node_Cmp
{
public:
	bool operator()(const N& a, const N&b)
	{
		return a.dist > b.dist;
	}
};

template<typename N>
class Edge
{
public:
	Edge(N*& inFrom, N*& inTo, typename N::Value_Type inLength) : from(inFrom), to(inTo), length(inLength), edge(Type::FORWARD)
	{
	}

	enum class Type
	{
		FORWARD,
		BACK,
		CROSS,
		TREE
	};

public:
	using Node_Type = typename N;
	N* from;
	N* to;
	enum class Type edge;
	typename N::Value_Type length;
};


template<typename T = Value<int32_t>>
class Graph
{
public:
	template<size_t sizex, size_t sizey>
	Graph(T(&inData)[sizex][sizey]) : size(sizex), data(new T*[size]), nodes(new Node<T>*[size])
	{
		static_assert(sizex == sizey, "Array must be of the square type");
		for (int x = 0; x < size; x++)
		{
			nodes[x] = new Node_Type(x);
			data[x] = new T[size];
		}

		for (int x = 0; x < size; x++)
		{
			for (int y = 0; y < size; y++)
			{
				data[x][y] = inData[x][y];
				if (data[x][y] != NO_EDGE)
					edges.emplace_back(new Edge_Type(nodes[x], nodes[y], data[x][y]));
			}
		}
	}

	Graph(T** (&&inData), size_t inSize) : size(inSize), data(inData), nodes(new Node<T>*[size])
	{
		for (int x = 0; x < size; x++)
			nodes[x] = new Node_Type(x);

		for (int x = 0; x < size; x++)
			for (int y = 0; y < size; y++)
				if (data[x][y] != NO_EDGE)
					edges.emplace_back(new Edge_Type(nodes[x], nodes[y], data[x][y]));

		inData = nullptr;
	}


	///*	 1	 2	 3	 4	 5	*/
	//	/*1*/{ z, 3, 3, z, z	 },
	//	/*2*/{ z, z, 4, 6, z },
	//	/*3*/{ z, z, z, -1, z },
	//	/*4*/{ 2, z, z, z, 1 },
	//	/*5*/{ -2, 3, 4, z, z }


	Graph<T>* addZeroNode(Node<T>*& outNode)
	{
		T** outData = new T*[size + 1];
		for (int x = 0; x < size + 1; x++)
			outData[x] = new T[size + 1];

		for (int x = 0; x < size; x++)
			for (int y = 0; y < size; y++)
				outData[x][y] = data[x][y];

		for (int x = 0; x < size + 1; x++)
		{
			outData[x][size] = NO_EDGE;
			outData[size][x] = 0;
		}
		outData[size][size] = NO_EDGE;

		Graph<T>* outGraph = new Graph<T>(move(outData), size + 1);
		outNode = outGraph->nodes[size];
		return outGraph;
	}

	Graph<T>* reverseEdges()
	{
		T** outData = new T*[size];
		for (int x = 0; x < size; x++)
			outData[x] = new T[size];

		for (int x = 0; x < size; x++)
		{
			for (int y = 0; y <= x; y++)
			{
				outData[x][y] = data[y][x];
				outData[y][x] = data[x][y];
			}
		}
		return new Graph<T>(move(outData), size);
	}

	~Graph()
	{
		for (int x = 0; x < size; x++)
			delete[] data[x];
		delete[] data;
		for (int x = 0; x < size; x++)
			delete nodes[x];
		delete[] nodes;
		for (auto& edge : edges)
			delete edge;
	}

public:
	using Node_Type = typename Node < T >;
	using Edge_Type = typename Edge < Node_Type >;
	static const T NO_EDGE;

	int size;
	T** data;
	Node_Type** nodes;
	vector<Edge_Type*> edges;
};

template <typename T>
const T Graph<T>::NO_EDGE = T::NONE;

template<typename G = Graph<>, typename N = G::Node_Type>
void initPaths(const G& graph, N*& start)
{
	for (int x = 0; x < graph->size; x++)
		graph->nodes[x]->reset();
	start->setDistance(0);
}

template<typename G = Graph<>, typename N = G::Node_Type>
void relax(N*& u, N*& v, G*& g)
{
	N::Value_Type newDist = u->dist + g->data[*u][*v];

	if (v->dist > newDist)
	{
		v->setDistance(newDist);
		v->setParent(*u);
	}
}

template<typename G = Graph<>, typename N = G::Node_Type>
bool dijkstra_relax(N& u, N& v, G*& g)
{
	N::Value_Type newDist = u.dist + g->data[u][v];

	if (v.dist > newDist)
	{
		v.setDistance(newDist);
		v.setParent(u);
		return true;
	}
	return false;
}

template<typename G = Graph<>, typename N = G::Node_Type>
void dijkstra(G*& graph, N*& start)
{
	initPaths(graph, start);
	priority_queue<N, vector<N>, Node_Cmp<N>> queue;
	for (int x = 0; x < graph->size; x++)
		queue.emplace(*graph->nodes[x]);
	while (!queue.empty())
	{
		N top = queue.top(); queue.pop();
		N* pTop = graph->nodes[top];
		if (pTop->state == N::State::CLOSED)
			continue;

		pTop->setDistance(top.dist);
		pTop->setParent(top.parent);
		pTop->close();

		for (int x = 0; x < graph->size; x++)
		{
			if (graph->data[top][x] != G::NO_EDGE)
			{
				N neigh = *graph->nodes[x];
				if (dijkstra_relax(top, neigh, graph))
					queue.emplace(neigh);
			}
		}
	}
}

template<typename G = Graph<>, typename N = G::Node_Type>
bool bellmanFord(G*& graph, N*& start)
{
	initPaths(graph, start);
	for (int x = 0; x < graph->size; x++)
		for (auto& edge : graph->edges)
			relax(edge->from, edge->to, graph);
	for (auto& edge : graph->edges)
		if (edge->to->dist > edge->from->dist + edge->length)
			return false;
	return true;
}

template<typename G = Graph<>, typename N = G::Node_Type>
bool johnson(G*& graph)
{
	N* zeroNode;
	G* zeroGraph = graph->addZeroNode(zeroNode);
	if (!bellmanFord(zeroGraph, zeroNode))
	{
		delete zeroGraph;
		return false;
	}

	auto delta = [zeroGraph](uint32_t inx) -> Value<int32_t>
	{
		static N::Value_Type zero(0);
		static const Value_Cmp<N::Value_Type> cmp;
		return min(zero, zeroGraph->nodes[inx]->dist, cmp);
	};

	for (auto& edge : graph->edges)
	{
		edge->length = edge->length + delta(*edge->from) - delta(*edge->to);
		graph->data[edge->from->index][edge->to->index] = edge->length;
	}

	for (int x = 0; x < graph->size; x++)
	{
		dijkstra(graph, graph->nodes[x]);
		for (int y = 0; y < graph->size; y++)
			zeroGraph->data[x][y] = graph->nodes[y]->dist - delta(x) + delta(y);
	}

	for (int x = 0; x < graph->size; x++)
		for (int y = 0; y < graph->size; y++)
			graph->data[x][y] = zeroGraph->data[x][y];
	delete zeroGraph;
	return true;
}

template<typename G = Graph<>, typename N = G::Node_Type>
void floydwarshall(G*& graph)
{
	Value_Cmp<N::Value_Type> cmp;
	int32_t size = graph->size;
	auto& d = graph->data;

	for (int k = 1; x <= graph->size; x++)
	{
		graph->nodes[x]->reset();
		graph->data[x][x] = 0;
	}

	for (int k = 1; x <= graph->size; x++)
		for (int i = 1; i <= graph->size; i++)
			for (int j = 1; j <= graph->size; j++)
				d[i][j] = min(d[i][j], d[i][k] + d[k][j], cmp);
}



template<typename G = Graph<>, typename N = G::Node_Type>
void dfs(G*& graph, N*& start)
{
	using S = N::State;
	stack<N*> stack;
	set<int32_t> nodes;
	int32_t time = 1;

	for (int x = 0; x < graph->size; x++)
		nodes.insert(x);

	while (!nodes.empty())
	{
		int32_t ind = *nodes.begin();
		stack.push(graph->nodes[ind]);
		while (!stack.empty())
		{
			N* cur = stack.top();
			nodes.erase(cur->index);

			switch (cur->state)
			{
			case S::OPEN:
				cur->state = S::CLOSED;
				cur->setCloseTime(time++); stack.pop();
				continue;
				break;
			case S::CLOSED:
				stack.pop();
				continue;
				break;
			default:
				cur->state = S::OPEN;
				cur->setOpenTime(time++);

				int32_t inx = cur->index;
				for (int x = graph->size - 1; x >= 0; x--)
					if (graph->data[inx][x] != G::NO_EDGE && graph->nodes[x]->state == S::FRESH)
						stack.push(graph->nodes[x]);
			}
		}
	}

	for (int x = 0; x < graph->size; x++)
		cout << (char)('1' + x) << ": " << setw(2) << graph->nodes[x]->oTime << " / " << setw(2) << graph->nodes[x]->cTime << endl;
	cout << endl;
};

int main(int argc, char* argv[])
{
	/*int32_t size = 5;
	int32_t** data = new int32_t*[size];
	for (int x = 0; x < size; x++)
		data[x] = new int32_t[size];

	srand((unsigned int)time(0));
	for (int x = 0; x < size; x++)
	{
		for (int y = 0; y <= x; y++)
		{
			if (x == y)
				data[x][y] = 0;
			else if (rand() % 3 == 1)
			{
				data[x][y] = ((rand() % (size - 1)) + 1);
				data[y][x] = 0;
			}
			else if (rand() % 3 == 2)
			{
				data[y][x] = ((rand() % (size - 1)) + 1);
				data[x][y] = 0;
			}
			else
			{
				data[y][x] = 0;
				data[x][y] = 0;
			}
		}
	}

	Graph<>* graph = new Graph<>(move(data), size);*/

	using V = Value < int32_t > ;
	const V z = Graph<V>::NO_EDGE;

	V data[][8] =
	{
		/*			 1	 2	 3	 4	 5	 6	 7	 8	 */
		/*1*/	{	 z,	 1,	 z,	 z,	 z,	 z,	 z,	 z,	 },
		/*2*/	{	 z,	 z,	 z,	 z,	 z,	 1,	 z,	 z,	 },
		/*3*/	{	 z,	 1,	 z,	 1,	 z,	 z,	 z,	 z,	 },
		/*4*/	{	 z,	 z,	 z,	 z,	 z,	 z,	 z,	 z,	 },
		/*5*/	{	 1,	 z,	 z,	 z,	 z,	 z,	 z,	 z,	 },
		/*6*/	{	 z,	 z,	 z,	 z,	 1,	 z,	 1,	 z,	 },
		/*7*/	{	 z,	 z,	 1,	 z,	 z,	 z,	 z,	 z,	 },
		/*8*/	{	 z,	 z,	 z,	 1,	 z,	 z,	 1,	 z,	 }
	};

	//V data[][5] =
	//{ 
	//			/*	 1	 2	 3	 4	 5	*/
	//	/*1*/	{	 z,	 3,	 3,	 z,	 z	 },
	//	/*2*/	{	 z,	 z,	 4,	 6,	 z	 },
	//	/*3*/	{	 z,	 z,	 z,	-1,	 z	 },
	//	/*4*/	{	 2,	 z,	 z,	 z,	 1	 },
	//	/*5*/	{	-2,	 3,	 4,	 z,	 z	 }
	//};

	Graph<V>* graph = new Graph<V>(data);
	//johnson(graph);
	dfs(graph, graph->nodes[0]);

	auto graph2 = graph->reverseEdges();
	dfs(graph2, graph2->nodes[0]);

	cout << "Finished" << endl;
	//int32_t test[][5] = 
	//{ 
	//			/*	1	2	3	4	5	*/
	//	/*1*/	{	0,	2,	3,	1,	0	},
	//	/*2*/	{	0,	0,	3,	2,	0	},
	//	/*3*/	{	-2, 0,	0,	2,	-1	},
	//	/*4*/	{	1,	0,	0,	0,	-2	},
	//	/*5*/	{	1,	3,	0,	0,	0	}
	//};

	//using T = Graph < int32_t > ;

	//int len[10] = { 0 };
	//for (int f = 0; f < 5; f++)
	//{
	//	T* graphA = new T(test);
	//	T* graphB = new T(test);

	//	dijkstra(graphA, graphA->nodes[f]);
	//	bellmanFord(graphB, graphB->nodes[f]);

	//	for (int x = 0; x < 5; x++)
	//	{
	//		if (x == f)
	//			continue;
	//		auto from = graphB->nodes[f];
	//		auto to = graphB->nodes[x];

	//		cout << "" << ((char)('1' + from->index)) << " -> " << ((char)('1' + to->index)) << ": ";
	//		cout << setw(2) << to->dist << endl;
	//		len[to->dist + 3]++;
	//	/*	if (to->dist != graphA->nodes[x]->dist)
	//		{
	//			cout << "Pruser" << endl;
	//			cout << "from: " << ((char)('1' + graphA->nodes[f]->index)) << " to: " << ((char)('1' + graphA->nodes[x]->index)) << ": ";
	//			cout << graphA->nodes[x]->dist << endl;
	//		}*/
	//	}
	//	//initPaths(*matrix, matrix->nodes[0]);
	//	//relax(matrix->nodes[1], matrix->nodes[0], *matrix);

	//		cout << endl;
	//	delete graphA;
	//	delete graphB;
	//}

	//cout << "distances: " << endl;
	//for (int x = 0; x < 8; x++)
	//{
	//	cout << setw(2) << (x - 3) << " x " << len[x] << endl;
 //	}

	system("pause");
	return 0;
}