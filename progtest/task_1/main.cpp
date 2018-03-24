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

int* parents;
class Edge* edges;
class Node* results;

class Edge
{
public:
	Edge() : from(0), to(0), price(0)
	{

	};

	int from;
	int to;
	int price;
};

class Node
{
public:
	set<int> friends;
};

struct cmpEdge
{
public:
	bool operator()(const Edge& lhs, const Edge& rhs)
	{
		return lhs.price < rhs.price;
	}
};

void makeSet(int x)
{
	parents[x] = x;
}


int findSet(int x)
{
	if (parents[x] != x)
		parents[x] = findSet(parents[x]);
	return parents[x];
}

void unionSet(int x, int y)
{
	int a = findSet(x);
	int b = findSet(y);

	if (a > b)
	{
		parents[a] = b;
	}
	else
	{
		parents[b] = a;
	}
}


void cd(const char * inFile, const char * outFile)
{
	ifstream input(inFile);
	int cities = 0;
	int pairs = 0;
	input >> cities >> pairs;

	parents = new int[cities];
	int* prices = new int[cities * cities];
	edges = new Edge[pairs];
	results = new Node[cities];

	for (int x = 0; x < cities * cities; x++)
		prices[x] = 0;

	for (int x = 0; x < cities; x++)
		makeSet(x);

	input.ignore(numeric_limits<streamsize>::max(), '\n');
	fstream::pos_type pos = input.tellg();
	input.ignore(numeric_limits<streamsize>::max(), '\n');

	int price;
	int edge = 0;
	int from, to;
	while (input >> from >> to >> price)
	{
		edges[edge].from = from;
		edges[edge].to = to;
		edges[edge].price = price;
		prices[from * cities + to] = price;
		prices[to * cities + from] = price;
		edge++;
	}

	long total = 0;
	input.clear();
	input.seekg(pos);
	char dummy;
	while (input >> from >> to >> dummy)
	{
		if (findSet(from) != findSet(to))
			unionSet(from, to);
		total += prices[from * cities + to];
		results[from].friends.insert(to);
		results[to].friends.insert(from);

		if (dummy != ',')
		{
			input.putback(dummy);
			break;
		}
	}
	input.close();

	sort(edges, edges + pairs, cmpEdge());

	for (int x = 0; x < pairs; x++)
	{
		int f = edges[x].from;
		int t = edges[x].to;
		if (findSet(f) != findSet(t))
		{
			total += prices[f * cities + t];
			results[f].friends.insert(t);
			results[t].friends.insert(f);
			unionSet(f, t);
		}
	}

	ofstream output(outFile);
	output << total << endl;
	for (int x = 0; x < cities; x++)
	{
		output << x;
		for (auto val : results[x].friends)
		{
			output << " " << val;
		} output << endl;
	}

	output.close();

	delete[] parents;
	delete[] results;
	delete[] edges;
	delete[] prices;
}


#ifndef __PROGTEST__
int main(int argc, char * argv[])
{
	cd("input.txt", "output.txt");
	return 0;
}
#endif /* __PROGTEST__ */
