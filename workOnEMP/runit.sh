#!/bin/bash
echo "start now"
g++ squareGenerateSolve.cpp  -lpthread -o execution

a=1
echo $a
while (($a<=50))
do
	echo "the $a time start"
	./execution
	let "a++"
done
