#include <deque>
#include <iostream>
#include <memory>
#include <ostream>
#include <raylib.h>
#include <raymath.h>
#include <string.h>
#include <string>
#include <vector>
#include "c_logger.h"

char log_buffer[1024];


#define RAYLIB_ENABLED 1
const char screenshot_path[] = "screenshots/take_2/";

Color green = {175, 200, 100, 255};
Color darkgreen = {45, 50, 20, 255};

#define GRID_CELL_SIZE 40

#define SCREEN_W 1920
#define SCREEN_H 1040
#define DEFAULT_NODES 0
#define USE_RANDOM_WEIGHT 0
#define ENABLE_SCREEN_CAPTURE 0 

//#define SCREEN_W 125
//#define SCREEN_H 105

#define SUR_D 20

const int8_t sur_d[SUR_D*2] = {0, -1, -1, 0, 1, 1, 1, 0, -1,  0, 0, 1, 1, 1, 0, -1, -1, -1};
const double weight_grid[] =  {0, 1, sqrt(2), 1, sqrt(2), 1, sqrt(2), 1, sqrt(2)};


Color COL_WHITE         = {255, 255, 255, 255};
Color COL_GRAY          = {100, 100, 100, 255};
Color COL_GREEN         = {15,  165, 10,  255};
Color COL_BRIGHT_GREEN  = {15,  255, 10,  255};
Color COL_DARK_GREEN    = {5,    85, 3,   255};
Color COL_RED           = {255, 100, 10,  255};
Color COL_BLUE          = {10,  100, 255, 255};
Color COL_MAGENTA       = {120, 10,  200, 255};
Color COL_BLACK         = {1,   1,   1,   255};

enum col_enum{
    _BLACK,
    _GREEN,
    _RED,
    _BLUE,
    _BGREEN,
    _DGREEN,
    _GRAY,
    _MAGENTA
};

Color col_array[] = {COL_BLACK, COL_GREEN, COL_RED, COL_BLUE, COL_BRIGHT_GREEN, COL_DARK_GREEN, COL_GRAY, COL_MAGENTA};

void xy_to_pos(int* pos, int x, int y, int columns){
    if (x < 0 || x >= columns || y < 0 || y >= columns){
        *pos = -2;
        return;
    }
    *pos = x + (y * columns);
    return;
}

void pos_to_xy(int* x, int* y, int pos, int columns){
    *x = pos % columns;
    *y = pos / columns;
    return;
}


class Edge{
    public:
        Edge(size_t index, size_t other_index, double weight)
            :m_index(index), o_index(other_index), m_weight(weight){} 
        Edge(const Edge& edge)
            : m_index(edge.m_index), m_weight(edge.m_weight){ } //std::cout<< "Edge copied." << '\n'; }
        Edge(const Edge&& other) noexcept
            : m_index(other.m_index), m_weight(other.m_weight){ } //std::cout<< "Edge Moved." << '\n'; }
 
        ~Edge(){}

        void set_index(size_t index){
            m_index = index;
        }

        const size_t get_index() const {
            return m_index;
        }

        size_t get_other_index(){ 
            return o_index;
        }

        const double get_weight() const {
            return m_weight;
        }

    private:
        size_t m_index;
        size_t o_index;
        double m_weight;
};

class Node{
    public:
        Node(){
        }
        ~Node(){}
        void set_index(size_t index){
            m_index = index;
        }

        const size_t get_index() const {
            return m_index;
        }

        void get_xy(int* x, int* y){
            *x = m_x;
            *y = m_y;
        }

        void set_xy(int x, int y){
            m_x = x;
            m_x = y;
        }

        const double get_weights(size_t other_node_index) const {
            //return m_edges.at(other_node_index).get_weight();
            //Not gonna work: other_node_index does not correspond to index in m_edges.
            //Need a translation unit to convert actual node indices to m_edges indices.
            return 0 * other_node_index;
        }

        const std::vector<Edge>& get_edges() const {
            return m_edges;
        }

        size_t get_edge_count() const{
            return get_edge_size();
        }

        size_t get_edge_size() const{
            return m_edges.size();
        }

        void connect(size_t index, double weight){
            m_edges.emplace_back(index, m_index, weight);
        }
        std::vector<Edge>& getm_edges(){
            return m_edges;
        }

    private:
        size_t m_index;
        std::vector<Edge> m_edges;
        int m_x;
        int m_y;
};

class Graph{
    public:
        Graph(size_t node_count) :
            m_node_count(node_count)
        {
            m_index = 0;
            m_total_edge_count = 0;
            m_nodes.reserve(m_node_count);
        }
        ~Graph(){}

        size_t create_node(int x, int y){
            size_t ix = m_nodes.size();
            m_nodes.push_back(std::make_unique<Node>());
            m_nodes.at(ix)->set_xy(x, y);
            m_nodes.at(ix)->set_index(ix);
            return ix;
        }

        void connect(size_t a, size_t b, double weight){
            if (a < m_nodes.size() && b < m_nodes.size()){
                m_nodes.at(a)->connect(b, weight);
                m_total_edge_count++;
            }
        }

        void connect_bi(size_t a, size_t b, double weight){
            if (a < m_nodes.size() && b < m_nodes.size()){
                m_nodes.at(a)->connect(b, weight);
                m_nodes.at(b)->connect(a, weight);
            }
        }
        size_t get_node_count() const {
            return m_nodes.size();
        }

        size_t get_edge_count() const {
            return m_total_edge_count;
        }

        void get_node_xy(size_t node_index, int* x, int* y) const {
            m_nodes.at(node_index)->get_xy(x, y);
        }

        double get_connection_weight(size_t node_index, size_t other_node_index) const {
            return m_nodes.at(node_index)->get_weights(other_node_index);
        }


        const std::vector<Edge>& get_node_edges(size_t node_index) const {
            return m_nodes.at(node_index)->get_edges();
        }


    private:
        size_t m_node_count;
        size_t m_index;
        std::vector<std::unique_ptr<Node>> m_nodes;
        size_t m_total_edge_count;
};

class Priority_Queue{
    public:
        Priority_Queue(){
            m_setup = false;
        }
        Priority_Queue(size_t item_count) : m_queue_size(item_count){
            m_arr.reserve(m_queue_size);
            m_setup = true;
        }
        ~Priority_Queue(){}

        bool is_empty() const {
            return m_arr.empty();
        }

        void setup(size_t item_count){
            m_queue_size = item_count;
            m_arr.reserve(m_queue_size);
            m_setup = true;
            sprintf(log_buffer, "Priority_Queue::%s: for %ld items.", __func__, m_queue_size);
            logger(log_buffer, 4);
        }

        void purge() {
            int _n;
            double _d;
            while(!is_empty()){
                pop(&_n, &_d);
            }
            return;
        }


        void insert(int node_index, double node_weight){
            Q_Node qn{node_index, node_weight};
            m_arr.emplace_back(qn);   //changed push_back to emplace_back
            shift_up(m_arr.size()-1);
        }
        void push(int node_index, double node_weight){
            insert(node_index, node_weight);
        }

        int pop(int* node_index, double* node_weight){
            int size = m_arr.size();
            int result = 0;
            if (size==0){
                result = -1;
                sprintf(log_buffer, "PQ::%s: called on empty queue.", __func__);
                logger(log_buffer, 4);
            }else{
                *node_index = m_arr[0].index;
                *node_weight = m_arr[0].weight;
                m_arr[0] = m_arr[size - 1];
                m_arr.pop_back();
                shift_down(0);
            }
            return result;
        }

        int get_min(int* node_index, double* node_weight){
            int result = 0;
            if (m_arr.empty()) {
                result =  -1;
            }else{
                *node_index = m_arr[0].index;
                *node_weight = m_arr[0].weight;
            }
            return result;
        }

    private:
        size_t m_queue_size;
        bool m_setup;

        int get_p(int i){ return (i-1) / 2;}
        int get_lc(int i){ return (i*2) + 1;}
        int get_rc(int i){ return (i*2) + 2;}
        struct Q_Node{
            int index;
            double weight;
        };
        std::vector<Q_Node> m_arr;

        void shift_up(int i){
            while (i>0 && m_arr[get_p(i)].weight > m_arr[i].weight){
                std::swap(m_arr[get_p(i)], m_arr[i]);
                i = get_p(i);
            }
        }

        void shift_down(size_t i){
            size_t n = m_arr.size();
            while (true){
                size_t min_index = i;
                size_t left = get_lc(i);
                size_t right = get_lc(i);
                if (left < n && m_arr[left].weight < m_arr[min_index].weight) 
                    min_index = left;
                if (right < n && m_arr[right].weight < m_arr[min_index].weight) 
                    min_index = right;

                if (min_index == i) {break; }
                std::swap(m_arr[i], m_arr[min_index]);
                i = min_index;
            }
        }

        void shift_down_recursive(int i){
            int min_index = i;
            size_t left = get_lc(i);
            size_t right = get_rc(i);
            if (left < m_arr.size() && m_arr[left].weight < m_arr[min_index].weight) 
                min_index = left;
            if (right < m_arr.size() && m_arr[right].weight < m_arr[min_index].weight) 
                min_index = right;
            if (i != min_index){
                std::swap(m_arr[i], m_arr[min_index]);
                shift_down(min_index);
            }
        }

};

#define MAX 99999.9

class A_Star{
    public:
        A_Star(Graph* graph) : 
            m_graph(graph) { }
        ~A_Star(){
            delete[] m_dist;
            delete[] m_prev;
            delete[] m_visited;
        }

        double heuristic(int node_a, int node_b){
            int xa, ya, xb, yb;
            int pos_a, pos_b;
            m_graph->get_node_xy(node_b, &xa, &ya);
            m_graph->get_node_xy(node_b, &xb, &yb);
            return abs(xa - xb) + abs(ya - yb);
        }
            

        int setup_path(size_t source_node, size_t destination_node){
            m_init();
            m_source_node = source_node;
            m_destination_node = destination_node;
            m_dist[m_source_node] = 0;
            const std::vector<Edge>& edge_list = m_graph->get_node_edges(m_source_node);
            for (const Edge& edge : edge_list){
                const size_t con_node = edge.get_index();
                const double con_weight = edge.get_weight();
                queue.push(con_node, con_weight);
                m_dist[con_node] = con_weight;
                m_prev[con_node] = m_source_node;
            }
            m_path_exist = -1;
            m_visited[source_node] = 1;
            return 1;
        }

        int step_path(int* colour_array){ 
            int return_value = 0;
            int visiting_node, r;
            double weight;
            if (!queue.is_empty()){
                r = queue.pop(&visiting_node, &weight);
                if (m_visited[visiting_node] == 0){
                    const std::vector<Edge>& v_edge_list = m_graph->get_node_edges(visiting_node);
                    for (const Edge& v_edge : v_edge_list){
                        const size_t v_node = v_edge.get_index();
                        double v_weight = v_edge.get_weight();
                        double pot_weight = m_dist[visiting_node]+v_weight;
                        if (colour_array[visiting_node] == _BLACK){
                            colour_array[visiting_node] = _DGREEN;
                        }
                        if (pot_weight < m_dist[v_node]){
                            m_dist[v_node] = pot_weight;
                            m_prev[v_node] = visiting_node;
                            double h_weight = pot_weight + heuristic(visiting_node, v_node);
                            queue.push(v_node, pot_weight + h_weight); //changed v_weight to pot_weight
                        }
                    }
                    m_visited[visiting_node] = 1;
                }
            }else{
                sprintf(log_buffer, "%s: Path found.", __func__);
                logger(log_buffer, 4);
                return_value = 1;
                if (m_dist[m_destination_node] < MAX){
                    m_path_exist = 1;
                }else{
                    m_path_exist = 0;
                }

            }
            return return_value;
        }

        // Returns 0 if marking in progress, 1 if completed.
        int mark_shortest_path(int* colour_array){ 
            int return_value = 1;
            if (m_path_exist != 1){
                return return_value;
            }
            if (m_dist[m_destination_node] < MAX){
                int* path = new int[m_node_count];
                size_t v = m_destination_node;
                int i = 0;
                while (v!=m_source_node){
                    //_graph.set_node_colour(v, blue);
                    if (colour_array[v] != _RED){
                        colour_array[v] = _BLUE;
                    }
                    path[i] = v;
                    v = m_prev[v];
                    i++;
                }
                path[i] = v;
                delete[] path;
            }
            return return_value;
        }


    private:
        Graph* m_graph;
        size_t m_node_count;
        Priority_Queue queue;
        double* m_dist;
        int* m_prev;
        bool* m_visited;
        size_t m_source_node;
        size_t m_destination_node;
        int m_path_exist;

        void m_init(){
            m_node_count = m_graph->get_node_count();
            m_dist = new double[m_node_count];
            m_prev = new int[m_node_count];
            m_visited = new bool[m_node_count];
            queue.setup(m_node_count);
            queue.purge();
            for (size_t i = 0; i < m_node_count; ++i) {
                m_dist[i] = MAX;
                m_prev[i] = -1;
                m_visited[i] = false;
            }
        }

};

class Dijkstra{
    public:
        Dijkstra(Graph* graph) : 
            m_graph(graph) { }
        ~Dijkstra(){
            delete[] m_dist;
            delete[] m_prev;
            delete[] m_visited;
        }

        int setup_path(size_t source_node, size_t destination_node){
            m_init();
            m_source_node = source_node;
            m_destination_node = destination_node;
            m_dist[m_source_node] = 0;
            const std::vector<Edge>& edge_list = m_graph->get_node_edges(m_source_node);
            for (const Edge& edge : edge_list){
                const size_t con_node = edge.get_index();
                const double con_weight = edge.get_weight();
                queue.push(con_node, con_weight);
                m_dist[con_node] = con_weight;
                m_prev[con_node] = m_source_node;
            }
            m_path_exist = -1;
            m_visited[source_node] = 1;
            return 1;
        }

        int step_path(int* colour_array){ 
            int return_value = 0;
            int visiting_node, r;
            double weight;
            if (!queue.is_empty()){
                r = queue.pop(&visiting_node, &weight);
                if (m_visited[visiting_node] == 0){
                    const std::vector<Edge>& v_edge_list = m_graph->get_node_edges(visiting_node);
                    for (const Edge& v_edge : v_edge_list){
                        const size_t v_node = v_edge.get_index();
                        double v_weight = v_edge.get_weight();
                        double pot_weight = m_dist[visiting_node]+v_weight;
                        if (colour_array[visiting_node] == _BLACK){
                            colour_array[visiting_node] = _DGREEN;
                        }
                        if (pot_weight < m_dist[v_node]){
                            m_dist[v_node] = pot_weight;
                            m_prev[v_node] = visiting_node;
                            queue.push(v_node, pot_weight); //changed v_weight to pot_weight
                        }
                    }
                    m_visited[visiting_node] = 1;
                }
            }else{
                sprintf(log_buffer, "%s: Path found.", __func__);
                logger(log_buffer, 4);
                return_value = 1;
                if (m_dist[m_destination_node] < MAX){
                    m_path_exist = 1;
                }else{
                    m_path_exist = 0;
                }

            }
            return return_value;
        }

        // Returns 0 if marking in progress, 1 if completed.
        int mark_shortest_path(int* colour_array){ 
            int return_value = 1;
            if (m_path_exist != 1){
                return return_value;
            }
            if (m_dist[m_destination_node] < MAX){
                int* path = new int[m_node_count];
                size_t v = m_destination_node;
                int i = 0;
                while (v!=m_source_node){
                    //_graph.set_node_colour(v, blue);
                    if (colour_array[v] != _RED){
                        colour_array[v] = _BLUE;
                    }
                    path[i] = v;
                    v = m_prev[v];
                    i++;
                }
                path[i] = v;
                delete[] path;
            }
            return return_value;
        }


    private:
        Graph* m_graph;
        size_t m_node_count;
        Priority_Queue queue;
        double* m_dist;
        int* m_prev;
        bool* m_visited;
        size_t m_source_node;
        size_t m_destination_node;
        int m_path_exist;

        void m_init(){
            m_node_count = m_graph->get_node_count();
            m_dist = new double[m_node_count];
            m_prev = new int[m_node_count];
            m_visited = new bool[m_node_count];
            queue.setup(m_node_count);
            queue.purge();
            for (size_t i = 0; i < m_node_count; ++i) {
                m_dist[i] = MAX;
                m_prev[i] = -1;
                m_visited[i] = false;
            }
        }

};

double select_weight(){
    double weight = 1.0;
    return weight;
}

class World{
    public:
        World(){}
        World(size_t screen_h, size_t screen_w, size_t cell_size) :
            m_screen_h(screen_h), m_screen_w(screen_w), m_cell_size(cell_size)
    {
        m_columns = screen_w / m_cell_size;
        m_rows    = screen_h / m_cell_size;
        m_cell_count = m_columns * m_rows;
        m_remainder_w = screen_w - (m_cell_size * m_columns);
        m_remainder_h = screen_h - (m_cell_size * m_rows);
        sprintf(log_buffer, "%s: height: %ld, width: %ld.", __func__, m_screen_h, m_screen_w);
        logger(log_buffer, 4);
        sprintf(log_buffer, "%s: cols: %ld, rows: %ld, cells: %ld.", __func__, m_columns, m_rows, m_cell_count);
        logger(log_buffer, 4);
        sprintf(log_buffer, "%s: rem_w: %ld, rem_h: %ld.", __func__, m_remainder_w, m_remainder_h);
        logger(log_buffer, 4);
        m_cell_colours = new int[m_cell_count];
        for (size_t c=0; c<m_cell_count; c++){
            m_cell_colours[c] = _BLACK;
        }
        m_square_green = -1;
        m_square_red = -1;


    }
        ~World(){
            delete[] m_cell_colours;
        }

        size_t get_cell_count(){
            return m_cell_count;
        }

        int get_cell_index_from_pos(int x, int y){
            size_t xp = (x - (m_remainder_w / 2)) / m_cell_size;
            size_t yp = (y - (m_remainder_h / 2)) / m_cell_size;
            int pos;
            xy_to_pos(&pos, xp, yp, m_columns);
            sprintf(log_buffer, "%s: x, y: (%d,%d), xyp: (%ld,%ld).", __func__, x, y, xp, yp);
            logger(log_buffer, 4);
            return pos;
        }

        void set_cell_colour(size_t pos, int colour){
            m_cell_colours[pos] = colour;
            if (colour == _RED){
                if (m_square_red != -1){
                    m_cell_colours[m_square_red] = _BLACK;
                }
                m_square_red = pos;
            }
            if (colour == _GREEN){
                if (m_square_green != -1){
                    m_cell_colours[m_square_green] = _BLACK;
                }
                m_square_green = pos;
            }
            if (colour == _BLACK){
                if (pos == m_square_green){
                    m_square_green = -1;
                }
                if (pos == m_square_red){
                    m_square_red = -1;
                }
            }
        }

        int get_cell_colour(size_t pos){
            return m_cell_colours[pos];
        }

        void draw_colour_cell(size_t x, size_t y, int colour){
            DrawRectangle(x+2, y+2, m_cell_size-2, m_cell_size-2, col_array[colour]);   
        }

        int construct_graph(){
            int ix;
            int x, y;
            for (size_t i=0; i<m_cell_count; i++){
                pos_to_xy(&x, &y, i, m_columns);
                ix = m_graph.create_node(x, y);
                if (ix != (int) i){
                    sprintf(log_buffer, "%s: node index error.", __func__ );
                    logger(log_buffer, 2);
                }
            }
            sprintf(log_buffer, "%s: Graph constructed with %ld nodes.", __func__, m_graph.get_node_count());
            logger(log_buffer, 4);
            return 1;

        }

        int connect_neighbours(){
            int x, y, n_node;
            for (size_t node=0; node<m_cell_count; node++){
                pos_to_xy(&x, &y, node, m_columns); // This node's position.
                for (int j=1; j<9; j++){
                    int sx = sur_d[j];
                    int sy = sur_d[j+9];
                    int ox = sx + x;
                    int oy = sy + y;

                    if (ox>=0 && oy>=0 && ox<(int)m_columns && oy<(int)m_rows ){
                        xy_to_pos(&n_node, ox, oy, m_columns); // Other node's postion.
                        if (n_node < 0){
                            sprintf(log_buffer, "%s: Node from i: %ld out of bounds: xy:(%d,%d), oxy:(%d,%d), sxy:(%d,%d).", __func__, node, x, y, ox, oy, sx, sy);
                            logger(log_buffer, 2);
                        }else{
                            if (get_cell_colour(n_node) != _MAGENTA){
                                double weight = weight_grid[j] * select_weight();
                                m_graph.connect(node, n_node, weight_grid[j] * weight);
                            }
                        }
                    }
                }
            }
            sprintf(log_buffer, "%s: Graph connected with %ld edges.", __func__, m_graph.get_edge_count());
            logger(log_buffer, 4);
            
            return 1;
        }

        int setup_astar(){
        //int setup_path(size_t source_node, size_t destination_node){
            if (m_square_green < 0 || m_square_red < 0){
                sprintf(log_buffer, "%s: No start or dist m_square selected.", __func__);
                logger(log_buffer, 2);
                std::cout<<log_buffer<<'\n';
            }else{
                m_astar.setup_path((size_t) m_square_green, (size_t) m_square_red);
            }

            //sprintf(log_buffer, "%s: completed.", __func__);
            //logger(log_buffer, 4);
            return 1;
        }

        int step_path_astar(){
            return m_astar.step_path(m_cell_colours);
        }

        int mark_path_astar(){
            return m_astar.mark_shortest_path(m_cell_colours);
        }



        int setup_dijkstra(){
        //int setup_path(size_t source_node, size_t destination_node){
            if (m_square_green < 0 || m_square_red < 0){
                sprintf(log_buffer, "%s: No start or dist m_square selected.", __func__);
                logger(log_buffer, 2);
                std::cout<<log_buffer<<'\n';
            }else{
                m_dijkstra.setup_path((size_t) m_square_green, (size_t) m_square_red);
            }

            //sprintf(log_buffer, "%s: completed.", __func__);
            //logger(log_buffer, 4);
            return 1;
        }

        int step_path_dijkstra(){
            return m_dijkstra.step_path(m_cell_colours);
        }

        int mark_path_dijkstra(){
            return m_dijkstra.mark_shortest_path(m_cell_colours);
        }

        void draw_cells(){
            int xp, yp;
            for (size_t i=0; i<m_cell_count; i++){
            //for (size_t i=0; i<4; i++){
                pos_to_xy(&xp, &yp, i, m_columns);
                size_t x = (xp * m_cell_size) + (m_remainder_w / 2);
                size_t y = (yp * m_cell_size) + (m_remainder_h / 2);
#if RAYLIB_ENABLED
                if (x < m_screen_w - m_remainder_w && y < m_screen_h - m_remainder_h){
                    DrawRectangleLines(x, y, m_cell_size, m_cell_size, COL_GRAY);
                    if (m_cell_colours[i] != _BLACK){
                        draw_colour_cell(x, y, m_cell_colours[i]);
                    }
                }
#endif
            }
        }


        int get_null(){
            return 0;
        }

            
    private:
        size_t m_screen_h;
        size_t m_screen_w;
        size_t m_cell_size;
        size_t m_columns;
        size_t m_remainder_w;
        size_t m_rows;
        size_t m_remainder_h;
        size_t m_cell_count;
        int* m_cell_colours;
        int m_square_green;
        int m_square_red;
        Graph m_graph = Graph(m_cell_size);
        Dijkstra m_dijkstra = Dijkstra(&m_graph);
        A_Star m_astar = A_Star(&m_graph);
};



void act_on_mouse(int mouse_button, World& world){
    int mouse_x = GetMouseX();
    int mouse_y = GetMouseY();
    int colour = _BLACK;
    if (mouse_button == 0){
        if (IsKeyDown(KEY_R)){
            colour = _RED;
        }
        if (IsKeyDown(KEY_B)){
            colour = _BLUE;
        }
        if (IsKeyDown(KEY_G)){
            colour = _GREEN;
        }
        if (IsKeyDown(KEY_M)){
            colour = _MAGENTA;
        }
    }else{
        colour = _BLACK;
    }
    int cell_index = world.get_cell_index_from_pos(mouse_x, mouse_y);
    world.set_cell_colour(cell_index, colour);
    sprintf(log_buffer, "%s: index: %d.", __func__, cell_index);
    logger(log_buffer, 4);
}


World world(SCREEN_H, SCREEN_W, GRID_CELL_SIZE);

void capture_screen(){
    static size_t capture_count = 0;
    char file_name[255];
    sprintf(file_name, "%sscreen_%ld.png", screenshot_path, capture_count);
    TakeScreenshot(file_name);
    sprintf(log_buffer, "%s: '%s'.", __func__, file_name);
    capture_count++;
    logger(log_buffer, 4);
}
    

int main(){

    int selector = 0;
    int (World::*world_actions[8])();
    world_actions[0] = &World::get_null;
    world_actions[1] = &World::construct_graph;
    world_actions[2] = &World::connect_neighbours;
    world_actions[3] = &World::setup_astar;
    world_actions[4] = &World::step_path_astar;
    world_actions[5] = &World::mark_path_astar;
    world_actions[6] = &World::get_null;

#if RAYLIB_ENABLED
    //InitWindow(GRID_CELL_SIZE*GRID_CELL_COUNT, GRID_CELL_SIZE*GRID_CELL_COUNT, "window");
    InitWindow(SCREEN_W, SCREEN_H, "Tiles");
    //SetTargetFPS(60);
    ToggleBorderlessWindowed();
    while (WindowShouldClose() == false){
        BeginDrawing();
        if (IsMouseButtonPressed(0)){
            act_on_mouse(0, world);
        }
        if (IsMouseButtonPressed(1)){
            act_on_mouse(1, world);
        }

        if(IsKeyPressed(KEY_ESCAPE) || IsKeyPressed(KEY_Q)){
            break;
        }
        if (IsKeyDown(KEY_LEFT_SHIFT) && IsKeyPressed(KEY_G)){
            //world.setup_dijkstra();
            if (selector == 0){
                selector = 1;
            }
        }
        if (DEFAULT_NODES){
            if (selector == 0){
                sprintf(log_buffer, "%s: system defined nodes.", __func__);
                logger(log_buffer, 4);
                world.set_cell_colour(0, _GREEN);
                world.set_cell_colour(world.get_cell_count()-1, _RED);
                world.setup_dijkstra();
                //selector++;
            }
        }
        selector += (world.*world_actions[selector])();

        if (selector < 7){
            ClearBackground(COL_BLACK);
            world.draw_cells();
        }
        EndDrawing();
        if (selector > 2 && selector <= 6){
            if (ENABLE_SCREEN_CAPTURE){
                capture_screen();
            }
        }
    }
    CloseWindow();
#endif //RL ENABLED

    sprintf(log_buffer, "Main is done.");
    logger(log_buffer, 4);

    return 0;
}

