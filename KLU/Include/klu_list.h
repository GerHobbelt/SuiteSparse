#ifndef _KLU_LIST_H_
#define _KLU_LIST_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

//#define int int

typedef struct node{
    struct node* next;
    struct node* prev;
    int value;
} node;

typedef struct list{
    node* head;
    node* tail;
    int length;
} list;

#define listsize sizeof(list)
#define nodesize sizeof(node)
#define BUBBLE 50
#define QUICK 51
#define FOUND 52
#define NOT_FOUND 53

void push_back(list*, int);
void push_front(list*, int);
void push_front_list(list*, list*);
void push_back_list(list*, list*);
int pop_front(list*);
int pop_back(list*);
void concat(list*, list*);
void concat_rev(list*, list*);
void sort_list(list*, int);
void bubblesort(list*);
void quicksort(list*, int);
void unlist(list*);
void killlist(list*);
void printlist(list*);
void printlistrev(list*);
int search(list*, const int);

#endif // _KLU_LIST_H_  
