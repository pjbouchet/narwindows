#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <string>
#include "parallel_hashmap/phmap.h"

int hash_lookup(std::int64_t cell){
  
  phmap::flat_hash_map<std::int64_t, int> geo = 
    {
    {306683210332, 1},
    {304553546288, 1},
    {308826245130, 1},
    {326214168375, 1}
    
    };
  
  return geo[cell];
}

#endif