
#include <iostream>
#include <chrono>
#include <string>
#include <queue>
#include<cmath>
#include <ilcplex/ilocplex.h>
#include <mysql_connection.h>
#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/prepared_statement.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>


using namespace std;
using namespace std::chrono;


// CONFIG ---------------------------------

string d = "8";
string intrsctProcess = "10";

// ----------------------------------------


class Equation {
public:
    int nVariables;
    std::vector<float> coefficients;
    char operand;
    float rightHandSide;

    Equation() {}

    Equation(int nVar, float coefs[], char op, float rhs) :
        nVariables(nVar),
        coefficients(coefs, coefs + nVar),
        operand(op),
        rightHandSide(rhs)
    {}

    Equation(const Equation& other) :
        nVariables(other.nVariables),
        coefficients(other.coefficients),
        operand(other.operand),
        rightHandSide(other.rightHandSide)
    {}

    Equation& operator=(const Equation& other) {
        if (this != &other) {
            nVariables = other.nVariables;
            coefficients = other.coefficients;
            operand = other.operand;
            rightHandSide = other.rightHandSide;
        }
        return *this;
    }

    void print() const {
        for (int i = 0; i < nVariables; i++) {
            std::cout << std::setw(8) << std::setprecision(3) << coefficients[i] << " X" << i + 1;
        }
        std::cout << std::setw(8) << operand << std::setw(8) << rightHandSide << std::endl;
    }
};


class Intersection {
public:
    int id;
    int f_i;
    int f_j;
    Equation eq;

    Intersection() {}

    Intersection(int iid, int fi, int fj, Equation equat)
        : id(iid), f_i(fi), f_j(fj), eq(equat) {}

    Intersection(const Intersection& other)
        : id(other.id), f_i(other.f_i), f_j(other.f_j), eq(other.eq) {}

    Intersection& operator=(const Intersection& other) {
        if (this != &other) {
            id = other.id;
            f_i = other.f_i;
            f_j = other.f_j;
            eq = other.eq;
        }
        return *this;
    }

    void set(int iid, int fi, int fj, Equation equat) {
        id = iid;
        f_i = fi;
        f_j = fj;
        eq = equat;
    }
};



class Node
{
public:
    int id=0;
    int f_i=0;
    int f_j=0;
    int intersection=0;
    Equation domain;
    int above=0;
    int below=0;
    int parent=0;
    int label=0;
    string proofs="";
    string hash="";
    float time=0;

    Node() {

    }

    Node(int fi, int fj, int intersect, Equation dom, int a, int b, int p, int l, string prf, string hsh, float t ) {
        f_i = fi;
        f_j = fj;
        intersection = intersect;
        domain = dom;
        above = a;
        below = b;
        parent = p;
        label = l;
        proofs = prf;
        hash = hsh;
        time = t;
    }

    void print() {
        cout << "id(" << id << ")  ";
        cout << "f_i(" << f_i << ")  ";
        cout << "f_j(" << f_j << ")  ";
        cout << "intersect(" << intersection << ")  ";
        cout << "above(" << above << ")  ";
        cout << "below(" << below << ")  ";
        cout << "parent(" << parent << ")  ";
        cout << "label(" << label << ")  ";
        cout << "time(" << time << ")"  << endl;
        cout << "Domain: ";
        domain.print();
    }
};



bool feasibilityChecking(Equation, Equation*, int);
void get_intersection(int, Intersection&);
void get_node(int, Node&);
int get_last_processed_intersection(void);
int get_num_of_intersections(void);
bool sql_connect(void);
void sql_disconnect(void);
string serialize(Equation);
std::vector<unsigned char> serializeBlob(Equation);
Equation deserialize(string);
int insert_node(Node);
bool update_above_below(int, int, int);
bool update_time(int, int, int, int);
int leftChildLabel(int);
int rightChildLabel(int);
int parentLabel(int);
string pathToRoot(int);
int pathToRootNumOfNodes(int);
void get_domain(int, Equation*);
int id2label(int);

std::unique_ptr<sql::ResultSet> sql_query(string);
sql::Driver* driver;
sql::Connection* con;
sql::Statement* stmt;


int main(int argc, char** argv)
{
    int numOfIntersect = 0;

    if (argc > 1) {
        numOfIntersect = atoi(argv[1]);
    }

    sql_connect();

    int lastIntersection = get_last_processed_intersection();

    // Check if ifmhtree table is empty --> insert first intersection as root
    if (!lastIntersection) {
        Intersection firstIntersection;
        get_intersection(1, firstIntersection);

        //Root
        float* coefs = new float[firstIntersection.eq.nVariables];
        for (int i = 0; i < firstIntersection.eq.nVariables; ++i) {
            coefs[i] = 1;
        }
        Equation D(firstIntersection.eq.nVariables, coefs, '>', 0);
        Node n(firstIntersection.f_i, firstIntersection.f_j, firstIntersection.id, D, 0, 0, 0, 1, "", "", 0);
        int id_of_root = insert_node(n);

        //Above
        firstIntersection.eq.operand = '>';
        Node a(firstIntersection.f_i, firstIntersection.f_j, firstIntersection.id, firstIntersection.eq, 0, 0, id_of_root, 2, "", "", 0);
        int id_of_a = insert_node(a);

        //Below
        firstIntersection.eq.operand = '<';
        Node b(firstIntersection.f_i, firstIntersection.f_j, firstIntersection.id, firstIntersection.eq, 0, 0, id_of_root, 3, "", "", 0);
        int id_of_b = insert_node(b);

        update_above_below(id_of_root, id_of_a, id_of_b);

        lastIntersection = 1;

        delete[] coefs;
    }

    

    int nEqs;
    Node currentNode;
    Intersection intersct;
    queue<int> nodesQueue;
    int numIntersections;
    if (numOfIntersect > 0) {
        numIntersections = numOfIntersect;
    }
    else {
        numIntersections = get_num_of_intersections();
    }
    cout << "Num of intersections: " << numIntersections << endl;
    cout << "Last processed intersection: " << lastIntersection << endl;

    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto durationPathToRoot = duration_cast<microseconds>(stop - start);
    auto durationFCandUpdate = duration_cast<microseconds>(stop - start);
    auto nonLeafProcessing = duration_cast<microseconds>(stop - start);


    for (int i = lastIntersection + 1; i < numIntersections; i++) {
        cout << "Intersection: " << i << endl;
        get_intersection(i, intersct);
        nodesQueue.push(1);
        while (nodesQueue.size() > 0) {
            get_node(nodesQueue.front(), currentNode);
            nodesQueue.pop();

            if (!currentNode.above && !currentNode.below) {
                start = high_resolution_clock::now();
                // Leaf
                nEqs = pathToRootNumOfNodes(currentNode.label);
                unique_ptr<Equation[]> eqs(new Equation[nEqs]);

                //getDomain
                std::unique_ptr<sql::ResultSet> res;
                res = sql_query("select domain from ifmhtree_"+ d + "d where find_in_set(label, '" + pathToRoot(id2label(currentNode.id)) + "')");
                int counter = 0;
                while (res->next()) {
                    if (!res->isNull(1)) eqs[counter] = deserialize(res->getString(1));
                    counter++;
                }
                stop = high_resolution_clock::now();
                durationPathToRoot = duration_cast<microseconds>(stop - start);
                start = high_resolution_clock::now();
                if (feasibilityChecking(intersct.eq, eqs.get(), nEqs)) {
                    //Above
                    intersct.eq.operand = '>';
                    Node a(intersct.f_i, intersct.f_j, intersct.id, intersct.eq, 0, 0, currentNode.id, leftChildLabel(currentNode.label), "", "", 0);
                    int id_of_a = insert_node(a);

                    //Below
                    intersct.eq.operand = '<';
                    Node b(intersct.f_i, intersct.f_j, intersct.id, intersct.eq, 0, 0, currentNode.id, rightChildLabel(currentNode.label), "", "", 0);
                    int id_of_b = insert_node(b);

                    update_above_below(currentNode.id, id_of_a, id_of_b);
                    
                    stop = high_resolution_clock::now();
                    durationFCandUpdate = duration_cast<microseconds>(stop - start);
                    update_time(currentNode.id, durationPathToRoot.count(), durationFCandUpdate.count(), nonLeafProcessing.count());
                    nonLeafProcessing = std::chrono::microseconds{ 0 };
                }
                else {
                    //cout << "No new partitions created" << endl;
                }
                eqs.reset(nullptr);
            }
            else {
                start = high_resolution_clock::now();
                // Non-leaf
                nEqs = pathToRootNumOfNodes(currentNode.label);
                unique_ptr<Equation[]> eqs(new Equation[nEqs]);

                //getDomain
                std::unique_ptr<sql::ResultSet> res;
                res = sql_query("select domain from ifmhtree_" + d + "d where find_in_set(label, '" + pathToRoot(id2label(currentNode.id)) + "')");
                int counter = 0;
                while (res->next()) {
                    if (!res->isNull(1)) eqs[counter] = deserialize(res->getString(1));
                    counter++;
                }

                if (feasibilityChecking(intersct.eq, eqs.get(), nEqs)) {
                    if (currentNode.above > 0) { nodesQueue.push(currentNode.above); }
                    if (currentNode.below > 0) { nodesQueue.push(currentNode.below); }
                }
                eqs.reset(nullptr);
                stop = high_resolution_clock::now();
                nonLeafProcessing = nonLeafProcessing + duration_cast<microseconds>(stop - start);
            }
        }
    }
}


bool sql_connect(void) {
    try {
        driver = get_driver_instance();
        con = driver->connect("tcp://127.0.0.1:3306", "root", "");
        con->setSchema("ifmh-tree");
        return true;
    }
    catch (sql::SQLException& e) {
        cout << "# ERR: SQLException in " << __FILE__;
        cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
        cout << "# ERR: " << e.what();
        cout << " (MySQL error code: " << e.getErrorCode();
        cout << ", SQLState: " << e.getSQLState() << " )" << endl;
        return false;
    }
}


void sql_disconnect(void) {
    
}

std::unique_ptr<sql::ResultSet> sql_query(string query) {
     try {
         stmt = con->createStatement();
         return std::unique_ptr<sql::ResultSet>(stmt->executeQuery(query));
     }
     catch (sql::SQLException& e) {
         cout << "# ERR: SQLException in " << __FILE__;
         cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
         cout << "# ERR: " << e.what();
         cout << " (MySQL error code: " << e.getErrorCode();
         cout << ", SQLState: " << e.getSQLState() << " )" << endl;
         return std::unique_ptr<sql::ResultSet>(nullptr);
     }
 }

void get_intersection(int id, Intersection& intersection) {
     int colsCount;
     std::unique_ptr<float[]> coefs;

     std::unique_ptr<sql::ResultSet> res = sql_query("select * from intersections_" + d + "d where id = '" + std::to_string(id) + "'");
     if (res->rowsCount() > 0) {
         res->next();
         colsCount = res->getMetaData()->getColumnCount(); //first 3 cols are not variables
         coefs.reset(new float[colsCount - 3]);
         for (int i = 4; i <= colsCount; i++) {
             coefs[i - 4] = stof(res->getString(i));
         }
         Equation eq(colsCount - 3, coefs.get(), '=', 0);
         intersection.set(stoi(res->getString(1)), stoi(res->getString(2)), stoi(res->getString(3)), eq);
     }
 }


 int get_last_processed_intersection(void) {
     //sql::ResultSet* res;
     std::unique_ptr<sql::ResultSet> res;
     res = sql_query("select max(intersection) from ifmhtree_" + d + "d");
     if (res->rowsCount() > 0) {
         res->next();
         if (!res->isNull(1))
             return stoi(res->getString(1));
         else
            return 0;
     }
     return 0;
 }



 int get_num_of_intersections(void) {
     //sql::ResultSet* res;
     std::unique_ptr<sql::ResultSet> res;
     res = sql_query("select count(*) from intersections_" + d + "d where id<='" + intrsctProcess+"'");
     if (res->rowsCount() > 0) {
         res->next();
         if (!res->isNull(1))
             return stoi(res->getString(1));
         else
             return 0;
     }
     return 0;
 }


 
 void get_node(int id, Node &node) {
     //static Node node;
     //sql::ResultSet* res;
     std::unique_ptr<sql::ResultSet> res;
     res = sql_query("select * from ifmhtree_" + d + "d where id = '" + std::to_string(id) + "'");
     if (res->rowsCount() > 0) {
         res->next();
         if (!res->isNull(1))  node.id = stoi(res->getString(1));
         if (!res->isNull(2))  node.f_i = stoi(res->getString(2));
         if (!res->isNull(3))  node.f_j = stoi(res->getString(3));
         if (!res->isNull(4))  node.intersection = stoi(res->getString(4));
         if (!res->isNull(5))  node.domain = deserialize(res->getString(5));
         if (!res->isNull(6))  node.above = stoi(res->getString(6));
         if (!res->isNull(7))  node.below = stoi(res->getString(7));
         if (!res->isNull(8))  node.parent = stoi(res->getString(8));
         if (!res->isNull(9))  node.label = stoi(res->getString(9));
     }
     //return node;
 }


int insert_node(Node node) {
     std::unique_ptr<sql::PreparedStatement> pstmt(con->prepareStatement("INSERT INTO ifmhtree_" + d + "d (f_i, f_j, intersection, domain, above, below, parent, label, proofs, hash, time1, time2, time3) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0, 0, 0)"));
     pstmt->setInt(1, node.f_i);
     pstmt->setInt(2, node.f_j);
     pstmt->setInt(3, node.intersection);
     pstmt->setString(4, serialize(node.domain));
     pstmt->setInt(5, node.above);
     pstmt->setInt(6, node.below);
     pstmt->setInt(7, node.parent);
     pstmt->setInt(8, node.label);
     pstmt->setString(9, node.proofs);
     pstmt->setString(10, node.hash);
     int rowsInserted = pstmt->executeUpdate();
     int id = 0;
     if (rowsInserted > 0) {
         std::unique_ptr<sql::Statement> stmt(con->createStatement());
         std::unique_ptr<sql::ResultSet> rs(stmt->executeQuery("SELECT LAST_INSERT_ID()"));
         if (rs->next()) {
             id = rs->getInt(1);
         }
     }
     return id;
 }

bool update_above_below(int id_root, int id_above, int id_below) {
     std::unique_ptr<sql::PreparedStatement> pstmt(con->prepareStatement("UPDATE ifmhtree_" + d + "d SET above = ? , below = ? WHERE id = ?"));
     pstmt->setInt(1, id_above);
     pstmt->setInt(2, id_below);
     pstmt->setInt(3, id_root);
     pstmt->executeUpdate();
     return true;
 }


 bool update_time(int id_root, int durationPathToRoot, int durationFCandUpdate, int nonLeafProcessing) {
     std::unique_ptr<sql::PreparedStatement> pstmt(con->prepareStatement("UPDATE ifmhtree_" + d + "d SET time1 = ?, time2 = ?, time3 = ?  WHERE id = ?"));
     pstmt->setInt(1, durationPathToRoot);
     pstmt->setInt(2, durationFCandUpdate);
     pstmt->setInt(3, nonLeafProcessing);
     pstmt->setInt(4, id_root);
     pstmt->executeUpdate();
     return true;
 }


 string serialize(Equation eq) {
     string serial = std::to_string(eq.nVariables) + "|";
     for (int i = 0; i < eq.nVariables; i++) {
         serial = serial + std::to_string(eq.coefficients[i]) + "|";
     }
     serial = serial + eq.operand + "|";
     serial = serial + std::to_string(eq.rightHandSide);
     return serial;
 }


 std::vector<unsigned char> serializeBlob(Equation eq) {
     std::vector<unsigned char> serial;
     serial.push_back(eq.nVariables);

     for (int i = 0; i < eq.nVariables; i++) {
         int coeff = eq.coefficients[i];
         serial.push_back((coeff >> 24) & 0xFF);
         serial.push_back((coeff >> 16) & 0xFF);
         serial.push_back((coeff >> 8) & 0xFF);
         serial.push_back(coeff & 0xFF);
     }

     serial.push_back(eq.operand);
     int rhs = eq.rightHandSide;
     serial.push_back((rhs >> 24) & 0xFF);
     serial.push_back((rhs >> 16) & 0xFF);
     serial.push_back((rhs >> 8) & 0xFF);
     serial.push_back(rhs & 0xFF);

     return serial;
 }

Equation deserialize(string eq) {
     std::unique_ptr<float[]> coefs;
     int nVariables = 0;
     int counter = 0;
     char operand;
     float rhs;

     const char delimiter = '|';
     int start = 0;

     string token;

     //Get the number of variables first --> otherwise error raise
     for (int i = 0; i < eq.length(); ++i) {
         if (eq[i] == delimiter) {
             token = eq.substr(start, i - start);
             break;
         }
     }

     nVariables = stoi(token);
     coefs = std::make_unique<float[]>(nVariables);

     for (int i = 0; i < eq.length(); ++i) {
         if (eq[i] == delimiter) {
             token = eq.substr(start, i - start);
             if (counter == 0) {
                 //nVariables is already taken --> Do nothing
             }
             else if (counter == nVariables + 1) {
                 operand = token[0];
             }
             else {
                 coefs[counter - 1] = stof(token);
             }
             counter++;
             start = i + 1;
         }
     }
     token = eq.substr(start, eq.length() - start);
     rhs = stof(token);

     Equation deserlial(nVariables, coefs.get(), operand, rhs);
     return deserlial;
 }



 int leftChildLabel(int id) {
    return pow(2, floor(log2(id)) + 1) + 2 * (id - pow(2, floor(log2(id))));
 } 


 int rightChildLabel(int id) {
     return leftChildLabel(id) +1;
 }


 int parentLabel(int id) {
     return (pow(2, floor(log2(id)) - 1) + ((id - pow(2, floor(log2(id)))) / 2));
 }
 
 
 string pathToRoot(int id) {
     int tmp = id;
     string path;
     
     while (tmp > 1) {
         path = path + std::to_string(tmp) + ",";
         tmp = parentLabel(tmp);
     }
     path = path + "1";
     return(path);
 }


 int pathToRootNumOfNodes(int id) {
     int tmp = id;
     int numOfNodes=0;

     while (tmp > 1) {
         numOfNodes++;
         tmp = parentLabel(tmp);
     }
     numOfNodes++;
     return(numOfNodes);
 }

 int id2label(int id) {
     //sql::ResultSet* res;
     std::unique_ptr<sql::ResultSet> res;
     res = sql_query("select label from ifmhtree_" + d + "d where id='" + std::to_string(id) + "'");
     if (res->next()) {
         if (!res->isNull(1)) {
             int result = stoi(res->getString(1));
             //delete res;
             return result;
         }
     }
     //delete res;
     return 0;
 }

void get_domain(int id, std::vector<Equation>& eqs) {
     int label = id2label(id);

     std::unique_ptr<sql::ResultSet> res;
     res = sql_query("select domain from ifmhtree_" + d + "d where find_in_set(label, '" + pathToRoot(label) + "')");

     // Get equations and fill eqs
     while (res->next()) {
         Equation eq = deserialize(res->getString(1));
         eq.print();
         eqs.push_back(eq);
     }
 }



bool feasibilityChecking(Equation target, Equation* eqs, int nEqs) {

    IloEnv   env;
    IloEnv   env2;
    try {

        //minimum ---------------------------------
        IloModel model(env);
        IloNumVarArray var(env);
        IloRangeArray con(env);

        IloExpr targ(env);

        for (int i = 0; i < target.nVariables; i++) {
            var.add(IloNumVar(env, 0, 1, ILOFLOAT));
            targ += target.coefficients[i] * var[i];
        }
        model.add(IloMinimize(env, targ));
        targ.end();
 

        for (int i = 0; i < nEqs; i++) {

            IloExpr condit(env);

            for (int j = 0; j < target.nVariables; j++) {
                condit += eqs[i].coefficients[j] * var[j];
            }

            if (eqs[i].operand == '>') {
                con.add(condit >= eqs[i].rightHandSide);
            }
            else if (eqs[i].operand == '<') {
                con.add(condit <= eqs[i].rightHandSide);
            }
            model.add(con);
            condit.end();
        }

        IloCplex cplex(model);
        cplex.setParam(cplex.PreInd, true);
        cplex.setParam(cplex.PreDual, true);

        cplex.setOut(env.getNullStream());

        if (!cplex.solve()) {
            throw(-1);
        }

        IloNumArray values(env);
        IloNumArray duals(env);

        cplex.getValues(values, var);
        cplex.getDuals(duals, con);





        //maximum ---------------------------------
        IloModel model2(env2);
        IloNumVarArray var2(env2);
        IloRangeArray con2(env2);

        IloExpr targ2(env2);

        for (int i = 0; i < target.nVariables; i++) {
            var2.add(IloNumVar(env2, 0, 1));
            targ2 += target.coefficients[i] * var2[i];
        }
        model2.add(IloMaximize(env2, targ2));
        targ2.end();

        for (int i = 0; i < nEqs; i++) {

            IloExpr condit2(env2);

            for (int j = 0; j < target.nVariables; j++) {
                condit2 += eqs[i].coefficients[j] * var2[j];
            }
            if (eqs[i].operand == '>') {
                con2.add(condit2 >= eqs[i].rightHandSide);
            }
            else if (eqs[i].operand == '<') {
                con2.add(condit2 <= eqs[i].rightHandSide);
            }
            model2.add(con2);
            condit2.end();
        }


        IloCplex cplex2(model2);

        cplex2.setOut(env2.getNullStream());

        if (!cplex2.solve()) {
            throw(-1);
        }

        IloNumArray values2(env2);
        IloNumArray duals2(env2);
        cplex2.getValues(values2, var2);
        cplex2.getDuals(duals2, con);


        if (cplex.getObjValue() < 0 && cplex2.getObjValue() > 0) {
            con.end();
            var.end();
            values.end();
            duals.end();
            con2.end();
            var2.end();
            values2.end();
            duals2.end();
            model.end();
            model2.end();
            env2.end();

            return true;
        }
        else {

            con.end();
            var.end();
            values.end();
            duals.end();
            con2.end();
            var2.end();
            values2.end();
            duals2.end();
            model.end();
            model2.end();
            env2.end();

            return false;
        }
    }
    catch (IloException& e) {
        env2.end();
        return false;
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        env2.end();
        return false;
        cerr << "Unknown exception caught" << endl;
    }

    return false;
}

