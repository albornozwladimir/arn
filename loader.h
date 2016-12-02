#ifndef _LOADER_H
#define _LOADER_H

#include <iostream>
#include <fstream>
#include <istream>
#include <math.h>
#include <string>
#include "constants.h"
#include "data.h"

void populate(const char *userdatadir,bool userdatalogic);
unsigned char getBase(std::string base);
unsigned char getBase1(std::string base);
int initStackValues(std::string fileName);
int initMiscloopValues(std::string fileName);
int initDangleValues(std::string fileName);
int initLoopValues(std::string fileName);
int initTstkhValues(std::string fileName);
int initTstkiValues(std::string fileName);
int initTloopValues(std::string fileName);
int initInt21Values(std::string fileName);
int initInt22Values(std::string fileName);
int initInt11Values(std::string fileName);

#endif
