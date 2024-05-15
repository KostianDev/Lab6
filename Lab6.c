#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <math.h>
#include <SDL2/SDL_ttf.h>
#include <stdbool.h>

#define n1 3
#define n2 1
#define n3 0
#define n4 3
#define N (10 + n3)
#define MAX_RAND 19
#define MIN_RAND 0
#define k (1.0 - n3 * 0.01 - n4 * 0.005 - 0.05)
#define WIDTH 1280
#define HEIGHT 720
#define SEED (1000 * n1 + 100 * n2 + 10 * n3 + n4)
#define ARROW_ANGLE  (M_PI / 5.0)
#define SIZE_MULT 3
#define SHIFT_ANGLE (M_PI / 18.0)

typedef struct passMatrix {
    double matrix[N][N];
} matrix;

typedef struct node {
    int key;
    int x;
    int y;
    int pos;
    bool used;
    struct node *next_node;
} l_list;

typedef struct Edge {
    int src;
    int dest;
    int weight;
    bool used;
    struct Edge *next_edge;
} e_list;

typedef struct Graph {
    bool finished;
    l_list *vertices;
    e_list *edges;
} graph;

bool KEYS[322];

l_list *l_list_init(int key, int x, int y, int pos);

l_list *addto_list(l_list *l_pointer, int key, int x, int y, int pos);

l_list *delfrom_start(l_list *l_pointer);

l_list *find_num(l_list *l_pointer, int key);

e_list *e_list_init(int src, int dest, int weight);

e_list *addto_edges(e_list *e_pointer, int src, int dest, int weight);

e_list *remove_edge(e_list *e_pointer);

e_list *swap(e_list *ptr1, e_list *ptr2);

int bubbleSort(e_list **head, int count);

void printList(FILE *fptr, e_list *n);

graph *graph_init(graph *graph1, int num_vertices, int num_edges);

matrix generateDirectedMatrix();

matrix generateUndirectedMatrix(matrix passMatrix);

matrix generateWeighMatrix(matrix undirectedMatrix);

void printMatrix(matrix passMatrix, FILE *fptr, const char[]);

void printGraph(graph *graph, FILE *fptr);

void clearScreen(SDL_Renderer *Renderer);

void drawCircle(SDL_Renderer *renderer, int32_t centreX, int32_t centreY, int32_t radius);

void drawBezierCurve(SDL_Renderer *renderer, int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4);

void drawArrowHead(SDL_Renderer *renderer, int endX, int endY, int gap, double angle);

void drawVertexNumber(SDL_Renderer *renderer, int number, int x, int y, int gap, SDL_Color color);

l_list *drawGraph(SDL_Renderer *renderer, SDL_Window *window, matrix graphMatrix, matrix weighMatrix,
                  l_list *list_ptr, int size, int isDir, int *r);

void drawDirConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                        int r, int dir, int width, int height, int gap2, int size);

void drawUndirConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                          int r, int width, int height, int size, matrix weightMatrix);

graph *primMST(graph *G, graph *G_T, int *isFirst, int *currV, int r,
               SDL_Renderer *renderer, SDL_Window *window, matrix weightMatrix, FILE *fptr);

int main(int argc, char *argv[]) {
    srand(SEED);

    FILE *fptr;

    TTF_Init();

    int num_vertices = N, init = 1, list_size = 0, flag = 1, isFirst = 1, vertex, r;

    fptr = fopen("Output.txt", "a");

    matrix directedMatrix = generateDirectedMatrix();
    matrix undirectedMatrix = generateUndirectedMatrix(directedMatrix);
    matrix weightMatrix = generateWeighMatrix(undirectedMatrix);

    printMatrix(undirectedMatrix, fptr, "undirected graph");
    printMatrix(weightMatrix, fptr, "weigh");
    fclose(fptr);

    l_list *G_ptr;

    e_list *edges_ptr;

    graph *G, *G_T;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (undirectedMatrix.matrix[i][j]){
                if (init) {
                    edges_ptr = e_list_init(i + 1, j + 1, (int) weightMatrix.matrix[i][j]);
                    init = 0;
                    list_size++;
                } else {
                    edges_ptr = addto_edges(edges_ptr, i + 1,
                                            j + 1, (int) weightMatrix.matrix[i][j]);
                    list_size++;
                }
            }
    printList(fptr, edges_ptr);

    bubbleSort(&edges_ptr, list_size);

    printList(fptr, edges_ptr);

    SDL_Init(SDL_INIT_EVERYTHING);

    SDL_Window *undirectedWindow = SDL_CreateWindow("Undirected Graph", SDL_WINDOWPOS_UNDEFINED,
                                                    SDL_WINDOWPOS_UNDEFINED, WIDTH,
                                                    HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer *undirectedRenderer = SDL_CreateRenderer(undirectedWindow, -1,
                                                          SDL_RENDERER_ACCELERATED);

    SDL_Event event;
    int quit = 0;

    while (!quit) {
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_WINDOWEVENT:
                    if (event.window.event == SDL_WINDOWEVENT_CLOSE) {
                        quit = 1;
                        break;
                    }
                case SDL_KEYDOWN:
                    if (event.key.keysym.sym < 322)
                        KEYS[event.key.keysym.sym] = true;
                    break;
                case SDL_KEYUP:
                    if (event.key.keysym.sym < 322)
                        KEYS[event.key.keysym.sym] = false;
                    break;
            }
        }

        if (flag) {
            clearScreen(undirectedRenderer);
            SDL_SetRenderDrawColor(undirectedRenderer, 255, 255, 255, 0);
            G_ptr = drawGraph(undirectedRenderer, undirectedWindow,
                              undirectedMatrix, weightMatrix, G_ptr, N, 0, &r);
            SDL_RenderPresent(undirectedRenderer);

            G = graph_init(G, N, list_size);
            G->vertices = G_ptr;
            G->edges = edges_ptr;

            G_T = graph_init(G_T, N, list_size);

            flag = 0;
        }

        if (KEYS[SDLK_n]) {
            G_T = primMST(G, G_T, &isFirst, &vertex, r, undirectedRenderer,
                          undirectedWindow, weightMatrix, fptr);
            KEYS[SDLK_n] = false;
        }

    }
    SDL_DestroyWindow(undirectedWindow);
    SDL_DestroyRenderer(undirectedRenderer);

    SDL_Quit();

    while (G_ptr != NULL)
        G_ptr = delfrom_start(G_ptr);

    while (edges_ptr != NULL)
        edges_ptr = remove_edge(edges_ptr);

    unlink("Output.txt");
    return 0;
}

l_list *l_list_init(int key, int x, int y, int pos) {
    l_list *l_pointer;
    l_pointer = malloc(sizeof(struct node));
    *l_pointer = (l_list) {
            .key = key,
            .x = x,
            .y = y,
            .pos = pos,
            .used = false,
            .next_node = NULL
    };
    return l_pointer;
}

l_list *addto_list(l_list *l_pointer, int key, int x, int y, int pos) {
    l_list *new_node;
    new_node = malloc(sizeof(struct node));
    new_node->key = key;
    new_node->x = x;
    new_node->y = y;
    new_node->pos = pos;
    new_node->used = false;
    new_node->next_node = l_pointer;
    return new_node;
}

l_list *delfrom_start(l_list *l_pointer) {
    l_list *node_ptr;
    node_ptr = l_pointer->next_node;
    free(l_pointer);
    return node_ptr;
}

l_list *find_num(l_list *l_pointer, int key) {
    l_list *this_node = l_pointer;

    while (this_node != NULL) {
        if (this_node->key == key) return this_node;
        else this_node = this_node->next_node;
    }
    return NULL;
}

e_list *e_list_init(int src, int dest, int weight) {
    e_list *e_pointer;
    e_pointer = malloc(sizeof(struct Edge));
    *e_pointer = (e_list) {
        .src = src,
        .dest = dest,
        .weight = weight,
        .used = false,
        .next_edge = NULL
    };
    return e_pointer;
}

e_list *addto_edges(e_list *e_pointer, int src, int dest, int weight) {
    e_list *new_edge;
    new_edge = malloc(sizeof(struct Edge));
    new_edge->src = src;
    new_edge->dest = dest;
    new_edge->weight = weight;
    new_edge->used = false;
    new_edge->next_edge = e_pointer;
    return new_edge;
}

e_list *remove_edge(e_list *e_pointer) {
    e_list *edge_ptr;
    edge_ptr = e_pointer->next_edge;
    free(e_pointer);
    return edge_ptr;
}

e_list *swap(e_list *ptr1, e_list *ptr2) {
    e_list *tmp = ptr2->next_edge;
    ptr2->next_edge = ptr1;
    ptr1->next_edge = tmp;
    return ptr2;
}

int bubbleSort(e_list **head, int count) {
    e_list **h;
    int i, j, swapped;

    for (i = 0; i <= count; i++) {

        h = head;
        swapped = 0;

        for (j = 0; j < count - i - 1; j++) {

            e_list *p1 = *h;
            e_list *p2 = p1->next_edge;

            if (p1->weight > p2->weight) {
                *h = swap(p1, p2);
                swapped = 1;
            }

            h = &(*h)->next_edge;
        }
        if (swapped == 0)
            break;
    }
}

void printList(FILE *fptr, e_list *n) {
    fptr = fopen("Output.txt", "a");
    fprintf(fptr, "\n\n");

    while (n != NULL) {
        fprintf(fptr, "%3.0d -> %3.0d     %d\n", n->src, n->dest, n->weight);
        n = n->next_edge;
    }
    fprintf(fptr, "\n");
    fclose(fptr);
}

graph *graph_init(graph *graph1, int num_vertices, int num_edges) {
    graph1 = malloc(sizeof(struct Edge) * num_edges + sizeof(struct node) * num_vertices);
    *graph1 = (graph) {
        .finished = false
    };
    return graph1;
}

matrix generateDirectedMatrix() {
    matrix passMatrix;
    double a;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            a = (rand() % (MAX_RAND + 1 - MIN_RAND) + MIN_RAND) / 10.0;
            a *= k;
            passMatrix.matrix[i][j] = a < 1.0 ? 0 : 1;
        }
    return passMatrix;
}

matrix generateUndirectedMatrix(matrix passMatrix) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (passMatrix.matrix[i][j] == 1) passMatrix.matrix[j][i] = 1;
    return passMatrix;
}

matrix generateWeighMatrix(matrix undirectedMatrix) {
    matrix B, C, D, H, Tr, W;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            B.matrix[i][j] = (rand() % (MAX_RAND + 1 - MIN_RAND) + MIN_RAND) / 100.0;
            C.matrix[i][j] = ceil(B.matrix[i][j] * 100 * undirectedMatrix.matrix[i][j]);
            if(C.matrix[i][j] == 0)
                D.matrix[i][j] = 0;
            else if(C.matrix[i][j] > 0)
                D.matrix[i][j] = 1;
        }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if(D.matrix[i][j] != D.matrix[j][i])
                H.matrix[i][j] = 1;
            else
                H.matrix[i][j] = 0;
        }

    for (int i = 0; i < N; ++i)
        for (int j = i; j < N; ++j) {
            Tr.matrix[j][i] = 0;
            Tr.matrix[i][j] = 1;
        }

    for (int i = 0; i < N; ++i)
        for (int j = i; j < N; ++j) {
            W.matrix[i][j] = (D.matrix[i][j] + H.matrix[i][j] * Tr.matrix[i][j]) * C.matrix[i][j];
            W.matrix[j][i] = W.matrix[i][j];
        }
    return W;
}

void printMatrix(matrix passMatrix, FILE *fptr, const char graphName[]) {
    fprintf(fptr, "\nHere is your %s matrix:\n\n", graphName);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(fptr, " %4.0lf", passMatrix.matrix[i][j]);
        }
        fprintf(fptr, "\n");
    }
    fprintf(fptr, "\n");
}

void printGraph(graph *graph, FILE *fptr) {
    int sumWeight = 0;

    fptr = fopen("Output.txt", "a");
    fprintf(fptr, "\nHere is your graph:\n\n");
    fprintf(fptr, "Vertices: ");

    l_list *tmpV = graph->vertices;
    e_list *tmpE = graph->edges;

    while (tmpV != NULL) {
        fprintf(fptr, " %d", tmpV->key);
        tmpV = tmpV->next_node;
    }

    fprintf(fptr, "\nEdges:\n");
    while (tmpE != NULL) {
        fprintf(fptr, "%3.0d -> %3.0d    %d\n", tmpE->src, tmpE->dest, tmpE->weight);
        sumWeight += tmpE->weight;
        tmpE = tmpE->next_edge;
    }
    fprintf(fptr, "\n\nSum of the edges' weights: %d", sumWeight);

    fclose(fptr);
}

void clearScreen(SDL_Renderer *Renderer) {
    SDL_SetRenderDrawColor(Renderer, 0, 0, 0, 255);
    SDL_RenderClear(Renderer);
}

void drawCircle(SDL_Renderer *renderer, int32_t centreX, int32_t centreY, int32_t radius) {
    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y) {
        SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

        if (error <= 0) {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0) {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }
}

void drawBezierCurve(SDL_Renderer *renderer, int x1, int x2, int x3, int x4, int y1, int y2, int y3, int y4) {
    double xu, yu, u;
    for (u = 0.0; u <= 1.0; u += 0.0005) {
        xu = pow(1 - u, 3) * x1 + 3 * u * pow(1 - u, 2) * x2 + 3 * pow(u, 2) * (1 - u) * x3
             + pow(u, 3) * x4;
        yu = pow(1 - u, 3) * y1 + 3 * u * pow(1 - u, 2) * y2 + 3 * pow(u, 2) * (1 - u) * y3
             + pow(u, 3) * y4;
        SDL_RenderDrawPoint(renderer, (int) xu, (int) yu);
    }
}

void drawArrowHead(SDL_Renderer *renderer, int endX, int endY, int gap, double angle) {
    double arrowSize = (double) gap / 3 / SIZE_MULT;

    double x1 = endX - arrowSize * cos(angle + ARROW_ANGLE);
    double y1 = endY - arrowSize * sin(angle + ARROW_ANGLE);
    double x2 = endX - arrowSize * cos(angle - ARROW_ANGLE);
    double y2 = endY - arrowSize * sin(angle - ARROW_ANGLE);

    SDL_Vertex vertex_1 = {{(float) endX, (float) endY},
                           {255,          255, 255, 255},
                           {1,            1}};
    SDL_Vertex vertex_2 = {{(float) x1, (float) y1},
                           {255,        255, 255, 255},
                           {1,          1}};
    SDL_Vertex vertex_3 = {{(float) x2, (float) y2},
                           {255,        255, 255, 255},
                           {1,          1}};

    SDL_Vertex vertices[] = {
            vertex_1,
            vertex_2,
            vertex_3
    };

    SDL_RenderGeometry(renderer, NULL, vertices, 3, NULL, 0);
}

void drawVertexNumber(SDL_Renderer *renderer, int number, int x, int y, int gap, SDL_Color color) {
    TTF_Font *font = TTF_OpenFont("arial.ttf", gap / 3);
    char numberString[4];

    snprintf(numberString, sizeof(numberString), "%d", number);

    SDL_Surface *surface = TTF_RenderText_Solid(font, numberString, color);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

    int texW, texH;
    SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);

    SDL_Rect rect = {x - texW / 2, y - texH / 2, texW, texH};
    SDL_RenderCopy(renderer, texture, NULL, &rect);

    SDL_FreeSurface(surface);
    SDL_DestroyTexture(texture);
    TTF_CloseFont(font);
}

l_list *drawGraph(SDL_Renderer *renderer, SDL_Window *window, matrix graphMatrix, matrix weighMatrix,
                  l_list *list_ptr, int size, int isDir, int *r) {
    int width, height, xTopCircles, xLowCircles, yCircles, turn = 1, mid, dir = 0;
    int gap, gap1, gap2, gap3, key = 1, pos = 1, degree = 0, outDegree = 0, inDegree = 0, flag = 1;

    if (size == 0) return NULL;

    SDL_GetWindowSize(window, &width, &height);
    SDL_Color white = {255, 255, 255, 255};

    if (size < 4) {
        if (size == 1) {
            *r = (height > width) ? width / 8 : height / 8;
            gap = *r * SIZE_MULT;

            drawCircle(renderer, width / 2, height / 2, *r);
            drawVertexNumber(renderer, key, width / 2, height / 2, gap, white);
            list_ptr = l_list_init(key, width / 2, height / 2, pos);
        } else {
            *r = (height > width) ? width / 12 : height / 12;

            gap = *r * SIZE_MULT;
            gap2 = (height - 2 * *r) / 2;
            gap3 = (width - 2 * gap - size * 2 * *r) / (size - 1);

            drawCircle(renderer, gap + *r, height / 2, *r);
            drawVertexNumber(renderer, key, gap + *r, height / 2, gap, white);
            list_ptr = l_list_init(key, gap + *r, height / 2, pos);
            key++;

            for (int i = 1; i < size; ++i) {
                drawCircle(renderer, gap + *r + i * (2 * *r + gap3),
                           height / 2, *r);
                drawVertexNumber(renderer, key, gap + *r + i * (2 * *r + gap3),
                                 height / 2, gap, white);
                list_ptr = addto_list(list_ptr, key,
                                      gap + *r + i * (2 * *r + gap3), height / 2, pos);
                key++;
            }
        }
    } else {

        mid = (abs(size - 5)) / 4;
        yCircles = 2 + mid * 2;
        xLowCircles = (size - yCircles) / 2;
        xTopCircles = size - xLowCircles - yCircles;

        *r = (height > width || (size >= 7 && size < 9)) ?
            width / (int) ((SIZE_MULT + 1) * (xTopCircles + 1) + SIZE_MULT) :
            height / ((SIZE_MULT + 1) * (yCircles / 2 + 1) + SIZE_MULT);

        gap = *r * SIZE_MULT;
        gap1 = (width - 2 * gap - (xTopCircles + 1) * 2 * *r) / xTopCircles;
        gap2 = (height - (yCircles / 2 + 1) * 2 * *r) / (yCircles / 2 + 2);
        gap3 = (width - 2 * gap - (xLowCircles + 1) * 2 * *r) / xLowCircles;

        drawCircle(renderer, gap + *r, gap2 + *r, *r);
        drawVertexNumber(renderer, key, gap + *r, gap2 + *r, gap, white);
        list_ptr = l_list_init(key, gap + *r, gap2 + *r, pos);
        key++;

        while (key <= size) {
            switch (turn % 4) {
                case 1:
                    for (int j = 1; j < xTopCircles; ++j) {
                        drawCircle(renderer, gap + *r + j * (2 * *r + gap1),
                                   gap2 + *r, *r);
                        drawVertexNumber(renderer, key, gap + *r + j * (2 * *r + gap1),
                                         gap2 + *r, gap, white);
                        list_ptr = addto_list(list_ptr, key,
                                              gap + *r + j * (2 * *r + gap1), gap2 + *r, pos);
                        key++;
                    }
                    break;
                case 2:
                    for (int j = 0; j < yCircles / 2; ++j) {
                        drawCircle(renderer, width - gap - *r,
                                   (gap2 + *r) + j * (gap2 + 2 * *r), *r);
                        drawVertexNumber(renderer, key, width - gap - *r,
                                         (gap2 + *r) + j * (gap2 + 2 * *r), gap, white);
                        list_ptr = addto_list(list_ptr, key,
                                              width - gap - *r, (gap2 + *r) + j * (gap2 + 2 * *r), pos);
                        key++;
                    }
                    break;
                case 3:
                    for (int j = 0; j < xLowCircles; ++j) {
                        drawCircle(renderer, width - gap - *r - j * (2 * *r + gap3),
                                   height - gap2 - *r, *r);
                        drawVertexNumber(renderer, key, width - gap - *r - j * (2 * *r + gap3),
                                         height - gap2 - *r, gap, white);
                        list_ptr = addto_list(list_ptr, key,
                                              width - gap - *r - j * (2 * *r + gap3),
                                              height - gap2 - *r, pos);
                        key++;
                    }
                    break;
                case 0:
                    for (int j = 0; j < yCircles / 2; ++j) {
                        drawCircle(renderer, gap + *r,
                                   height - (*r + gap2) - j * (2 * *r + gap2), *r);
                        drawVertexNumber(renderer, key, gap + *r,
                                         height - (*r + gap2) - j * (2 * *r + gap2), gap, white);
                        list_ptr = addto_list(list_ptr, key,
                                              gap + *r, height - (*r + gap2) - j * (2 * *r + gap2), pos);
                        key++;
                    }
                    break;
            }
            turn++;
            pos++;
        }
    }

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            if (graphMatrix.matrix[i][j] == 1) {
                l_list *node1 = find_num(list_ptr, i + 1);
                l_list *node2 = find_num(list_ptr, j + 1);
                if (i == j) {
                    if (node1->pos >= 3) flag = -1;
                    drawCircle(renderer,
                               (int) (node1->x - flag * (*r + (double) *r / 2 - flag)
                                                 * cos(M_PI / 4)),
                               (int) (node1->y - flag * (*r + (double) *r / 2 - flag)
                                                 * sin(M_PI / 4)),
                               *r / 2);
                    drawArrowHead(renderer, (int) (node1->x - flag * *r * cos(M_PI / 4)),
                                  (int) (node1->y - flag * *r * sin(M_PI / 4)), gap,
                                  3 * M_PI / 4 - flag * M_PI / 15);
                } else {
                    if(!isDir){
                        graphMatrix.matrix[j][i] = 0;
                        drawUndirConnections(renderer, node1, node2, *r, width, height, size, weighMatrix);
                    } else {
                        dir = 1;
                        drawDirConnections(renderer, node1, node2, *r, dir, width, height, gap2, size);
                    }
                    if (graphMatrix.matrix[j][i] == 1 && isDir != 0) {
                        graphMatrix.matrix[j][i] = 0;
                        dir = -1;
                        drawDirConnections(renderer, node2, node1, *r, dir, width, height, gap2, size);
                    }
                }
            }
    return list_ptr;
}

void drawDirConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                        int r, int dir, int width, int height, int gap2, int size) {
    int startX = node1->x, startY = node1->y, endX = node2->x, endY = node2->y;
    int mid1X, mid1Y, mid2X, mid2Y;
    int midX = (node1->x + node2->x) / 2;
    int midY = (node1->y + node2->y) / 2;
    int dirX = (node1->x - width / 2) <= 0 ? -1 : 1;
    int dirY = (node1->y - height / 2) <= 0 ? -1 : 1;
    int gap = SIZE_MULT * r;
    double angle1, angle2;
    if ((abs(node1->key - node2->key) == 1) ||
        (((node1->key == 1 || node2->key == 1) && (node1->key == size || node2->key == size)) && size >= 4)) {
        angle1 = atan2(endY - startY, endX - startX);
        SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1 + SHIFT_ANGLE)),
                           startY + (int) ((double) r * sin(angle1 + SHIFT_ANGLE)),
                           endX - (int) ((double) r * cos(angle1 - SHIFT_ANGLE)),
                           endY - (int) ((double) r * sin(angle1 - SHIFT_ANGLE)));
        drawArrowHead(renderer, endX - (int) ((double) r * cos(angle1 - SHIFT_ANGLE)),
                      endY - (int) ((double) r * sin(angle1 - SHIFT_ANGLE)), gap, angle1);
    } else {
        if (startX == endX && (startX == gap + r || startX == width - gap - r)) {
            mid1X = (int) ((double) (startX + midX) / 2 + gap * dir * dirX);
            mid1Y = (startY + midY) / 2;
            mid2X = (int) ((double) (midX + endX) / 2 + gap * dir * dirX);
            mid2Y = (midY + endY) / 2;
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                            mid2X, endX - (int) ((double) r * cos(angle2)),
                            startY + (int) ((double) r * sin(angle1)), mid1Y,
                            mid2Y, endY - (int) ((double) r * sin(angle2)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle2)),
                          endY - (int) ((double) r * sin(angle2)), gap, angle2);
        } else if (startY == endY && (startY == gap2 + r || startY == height - gap2 - r)) {
            mid1X = (startX + midX) / 2;
            mid1Y = (int) ((double) (startY + midY) / 2 + gap * dir * dirY);
            mid2X = (midX + endX) / 2;
            mid2Y = (int) ((double) (midY + endY) / 2 + gap * dir * dirY);
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                            mid2X, endX - (int) ((double) r * cos(angle2)),
                            startY + (int) ((double) r * sin(angle1)), mid1Y,
                            mid2Y, endY - (int) ((double) r * sin(angle2)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle2)),
                          endY - (int) ((double) r * sin(angle2)), gap, angle2);
        } else {
            angle1 = atan2(endY - startY, endX - startX);
            SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1 + SHIFT_ANGLE)),
                               startY + (int) ((double) r * sin(angle1 + SHIFT_ANGLE)),
                               endX - (int) ((double) r * cos(angle1 - SHIFT_ANGLE)),
                               endY - (int) ((double) r * sin(angle1 - SHIFT_ANGLE)));
            drawArrowHead(renderer, endX - (int) ((double) r * cos(angle1 - SHIFT_ANGLE)),
                          endY - (int) ((double) r * sin(angle1 - SHIFT_ANGLE)), gap, angle1);
        }
    }
}

void drawUndirConnections(SDL_Renderer *renderer, l_list *node1, l_list *node2,
                          int r, int width, int height, int size, matrix weightMatrix) {
    int startX = node1->x, startY = node1->y, endX = node2->x, endY = node2->y;
    int mid1X, mid1Y, mid2X, mid2Y;
    int midX = (node1->x + node2->x) / 2;
    int midY = (node1->y + node2->y) / 2;
    int dirX = (node1->x - width / 2) <= 0 ? -1 : 1;
    int dirY = (node1->y - height / 2) <= 0 ? -1 : 1;
    int gap = SIZE_MULT * r, yCircles, mid;

    if (size >= 4) {
        mid = (abs(size - 5)) / 4;
        yCircles = 2 + mid * 2;
    }

    int gap2 = (size < 4) ? (height - 2 * r) / 2 :
               (height - (yCircles / 2 + 1) * 2 * r) / (yCircles / 2 + 2);
    double angle1, angle2;

    SDL_Color color = {255, 255, 0, 255};
    TTF_Font *font = TTF_OpenFont("arial.ttf", gap / 7);
    char numberString[4];
    int weighX, weighY;

    snprintf(numberString, sizeof(numberString), "%d",
             (int) weightMatrix.matrix[node1->key - 1][node2->key - 1]);

    SDL_Surface *surface = TTF_RenderText_Solid(font, numberString, color);
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

    int texW, texH;
    SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);

    if ((abs(node1->key - node2->key) == 1) ||
        ((node1->key == 1 || node2->key == 1) && (node1->key == size || node2->key == size))) {
        angle1 = atan2(endY - startY, endX - startX);
        SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1)),
                           startY + (int) ((double) r * sin(angle1)),
                           endX - (int) ((double) r * cos(angle1)),
                           endY - (int) ((double) r * sin(angle1)));
        weighX =  startX + (int) ((double) r * 2 * cos(angle1));
        weighY = startY + (int) ((double) r * 2 * sin(angle1));
        SDL_Rect rect = {weighX - texW / 2, weighY - texH / 2, texW, texH};
        SDL_RenderCopy(renderer, texture, NULL, &rect);
    } else {
        if (startX == endX && (startX == gap + r || startX == width - gap - r)) {
            mid1X = (int) ((double) (startX + midX) / 2 + gap * dirX);
            mid1Y = (startY + midY) / 2;
            mid2X = (int) ((double) (midX + endX) / 2 + gap * dirX);
            mid2Y = (midY + endY) / 2;
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                        mid2X, endX - (int) ((double) r * cos(angle2)),
                        startY + (int) ((double) r * sin(angle1)), mid1Y,
                        mid2Y, endY - (int) ((double) r * sin(angle2)));
            weighX = (startX + (int) ((double) r * cos(angle1)) + mid1X) / 2;
            weighY = (startY + (int) ((double) r * sin(angle1)) + mid1Y) / 2;
            SDL_Rect rect = {weighX - texW / 2, weighY - texH / 2, texW, texH};
            SDL_RenderCopy(renderer, texture, NULL, &rect);
        } else if (startY == endY && (startY == gap2 + r || startY == height - gap2 - r)) {
            mid1X = (startX + midX) / 2;
            mid1Y = (int) ((double) (startY + midY) / 2 + gap * dirY);
            mid2X = (midX + endX) / 2;
            mid2Y = (int) ((double) (midY + endY) / 2 + gap * dirY);
            angle1 = atan2(mid1Y - startY, mid1X - startX);
            angle2 = atan2(endY - mid2Y, endX - mid2X);
            drawBezierCurve(renderer, startX + (int) ((double) r * cos(angle1)), mid1X,
                        mid2X, endX - (int) ((double) r * cos(angle2)),
                        startY + (int) ((double) r * sin(angle1)), mid1Y,
                        mid2Y, endY - (int) ((double) r * sin(angle2)));
            weighX = (startX + (int) ((double) r * cos(angle1)) + mid1X) / 2;
            weighY = (startY + (int) ((double) r * sin(angle1)) + mid1Y) / 2;
            SDL_Rect rect = {weighX - texW / 2, weighY - texH / 2, texW, texH};
            SDL_RenderCopy(renderer, texture, NULL, &rect);
        } else {
            angle1 = atan2(endY - startY, endX - startX);
            SDL_RenderDrawLine(renderer, startX + (int) ((double) r * cos(angle1)),
                               startY + (int) ((double) r * sin(angle1)),
                               endX - (int) ((double) r * cos(angle1)),
                               endY - (int) ((double) r * sin(angle1)));
            weighX = startX + (int) ((double) r * 2.5 * cos(angle1));
            weighY = startY + (int) ((double) r * 2.5 * sin(angle1));
            SDL_Rect rect = {weighX - texW / 2, weighY - texH / 2, texW, texH};
            SDL_RenderCopy(renderer, texture, NULL, &rect);
        }
    }

    SDL_FreeSurface(surface);
    SDL_DestroyTexture(texture);
    TTF_CloseFont(font);
}

graph *primMST(graph *G, graph *G_T, int *isFirst, int *currV, int r,
               SDL_Renderer *renderer, SDL_Window *window, matrix weightMatrix, FILE *fptr) {

    if (G_T->finished)
        return G_T;

    int gap = r * SIZE_MULT, width, height;

    l_list *tmpDest, *tmpSrc;
    e_list *tmpE;

    SDL_Color green = {0, 255, 0, 255};
    SDL_GetWindowSize(window, &width, &height);

    if (*isFirst == 1) {
        int a = rand() % (10 + 1 - 1) + 1;
        l_list *tmp = find_num(G->vertices, a);
        G_T->vertices = l_list_init(a, tmp->x, tmp->y, tmp->pos);
        tmp->used = true;
        *currV = a;

        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
        drawCircle(renderer, G_T->vertices->x, G_T->vertices->y, r);
        drawVertexNumber(renderer, G_T->vertices->key, G_T->vertices->x, G_T->vertices->y, gap, green);
        SDL_RenderPresent(renderer);

        *isFirst = 2;
        return G_T;
    }

    tmpE = G->edges;

    while (tmpE != NULL) {
        tmpDest = find_num(G->vertices, tmpE->dest);
        if (!tmpE->used && tmpE->src == *currV && !tmpDest->used && tmpE->src != tmpE->dest) {
            if (*isFirst == 2) {
                G_T->edges = e_list_init(*currV, tmpDest->key, tmpE->weight);
            } else {
                G_T->edges = addto_edges(G_T->edges, *currV, tmpDest->key, tmpE->weight);
            }
            G_T->vertices = addto_list(G_T->vertices, tmpDest->key, tmpDest->x, tmpDest->y, tmpDest->pos);
            tmpE->used = true;
            tmpDest->used = true;

            SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
            drawUndirConnections(renderer, find_num(G->vertices, *currV), tmpDest,
                                 r, width, height, N, weightMatrix);
            drawCircle(renderer, tmpDest->x, tmpDest->y, r);
            drawVertexNumber(renderer, tmpDest->key, tmpDest->x, tmpDest->y, gap, green);
            SDL_RenderPresent(renderer);

            *currV = tmpDest->key;
            *isFirst = 0;
            return G_T;
        } else tmpE = tmpE->next_edge;
    }

    tmpE = G->edges;

    while (tmpE != NULL) {
        tmpSrc = find_num(G->vertices, tmpE->src);
        tmpDest = find_num(G->vertices, tmpE->dest);
        if (!tmpE->used && !tmpDest->used && tmpE->src != tmpE->dest) {
            G_T->edges = addto_edges(G_T->edges, tmpSrc->key, tmpDest->key, tmpE->weight);
            G_T->vertices = addto_list(G_T->vertices, tmpDest->key, tmpDest->x, tmpDest->y, tmpDest->pos);
            tmpE->used = true;
            tmpDest->used = true;

            SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
            drawUndirConnections(renderer, tmpSrc, tmpDest, r, width, height, N, weightMatrix);
            drawCircle(renderer, tmpDest->x, tmpDest->y, r);
            drawVertexNumber(renderer, tmpDest->key, tmpDest->x, tmpDest->y, gap, green);
            SDL_RenderPresent(renderer);

            *currV = tmpDest->key;
            return G_T;
        } else tmpE = tmpE->next_edge;
    }

    G_T->finished = true;
    printGraph(G_T, fptr);
    return G_T;
}