#ifndef __PROGTEST__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>
#include <list>
#include <stdint.h>

using namespace std;
#endif /* __PROGTEST__ */

class Node
{
public:
	enum class State
	{
		FRESH,
		OPEN,
		CLOSED
	};

	Node() : state(State::FRESH), comp(0)
	{

	}

	State state;
	int comp;
	vector<int> friends;
};


Node* graph = nullptr;
Node* inverse = nullptr;
int cnt;
list<int> closed;
int comp;

void DFS_INVERSE_REC(int x)
{
	list<int> stack;
	stack.push_back(x);

	while (!stack.empty())
	{
		Node& cur = inverse[stack.back()];
		cur.comp = comp;
		if (cur.state == Node::State::OPEN)
		{
			cur.state = Node::State::CLOSED;
			stack.pop_back();
			continue;
		}
		else if (cur.state == Node::State::CLOSED)
		{
			stack.pop_back();
			continue;
		}
		cur.state = Node::State::OPEN;
		for (int neigh : cur.friends)
		{
			if (inverse[neigh].state == Node::State::FRESH)
				stack.push_back(neigh);
		}
	}
}

void DFS_GRAPH_REC(int x)
{
	list<int> stack;
	stack.push_back(x);

	while (!stack.empty())
	{
		Node& cur = graph[stack.back()];
		if (cur.state == Node::State::OPEN)
		{
			cur.state = Node::State::CLOSED;
			closed.push_front(stack.back());
			stack.pop_back();
			continue;
		}
		else if (cur.state == Node::State::CLOSED)
		{
			stack.pop_back();
			continue;
		}
		cur.state = Node::State::OPEN;
		for (int neigh : cur.friends)
		{
			if (graph[neigh].state == Node::State::FRESH)
				stack.push_back(neigh);
		}
	}
}

void DFS_GRAPH()
{
	for (int x = 0; x < cnt; x++)
	{
		if (graph[x].state == Node::State::FRESH)
			DFS_GRAPH_REC(x);
	}
}

void fit(const char * inFile, const char * outFile)
{
	closed.clear();
	ifstream input(inFile);
	cnt = 0;
	input >> cnt;

	int from, to;
	graph = new Node[cnt];
	inverse = new Node[cnt];

	while (input >> from >> to)
	{
		graph[from].friends.push_back(to);
		inverse[to].friends.push_back(from);
	}
	input.close();

	DFS_GRAPH();


	comp = 1;
	for (auto start : closed)
	{
		if (inverse[start].state == Node::State::FRESH)
		{
			DFS_INVERSE_REC(start);
			comp++;
		}
	}

	ofstream output(outFile);
	for (int x = 0; x < cnt; x++)
	{
		int cmp1 = inverse[x].comp;
		for (auto fr : graph[x].friends)
		{
			if (cmp1 != inverse[fr].comp)
			{
				output << x << endl;
				break;
			}
		}
	}
	output.close();

	delete[] graph;
	delete[] inverse;
}


#ifndef __PROGTEST__
int main(int argc, char * argv[])
{
	fit("input_a.txt", "output_a.txt");
	fit("input_b.txt", "output_b.txt");
	// system("pause");
	return 0;
}
#endif /* __PROGTEST__ */