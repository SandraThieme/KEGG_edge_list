#!/bin/bash

for i in `ls *.xml`; do echo $i; 

kegfile=$i

sed 's/</ /g' <"$kegfile" >test1.txt
sed 's/>/ /g' <test1.txt >test2.txt
sed 's_/_ _g' <test2.txt >test1.txt
kegfile="test1.txt"

while read -r prefix content
do
    case "$prefix" in
	entry) echo -e "$content";;
    esac
done < "$kegfile" >test3.txt

sed 's/ /;/g' <test3.txt >entry.txt
sed 's/id=/ /' <entry.txt >test3.txt
sed 's/name=/ /' <test3.txt >entry.txt
sed 's/type=/ /' <entry.txt >test3.txt
sed 's/reaction=/ /' <test3.txt >test2.txt
sed 's/ //g' <test2.txt >entry$i.txt

while read -r prefix content
do
    case "$prefix" in
	reaction) col1="$content";;
	substrate) col2="$content";;
	product) echo -e "$col1$col2$content";;
    esac
done < "$kegfile" >reaction.txt

sed 's/name=/ /g' <reaction.txt >test2.txt
sed 's/id=/ /g' <test2.txt >test3.txt
sed 's/type=/ /g' <test3.txt >test2.txt
sed 's/ /;/g' <test2.txt >test3.txt
sed 's/;;/;/g' <test3.txt >reaction$i.txt

rm test1.txt test2.txt test3.txt

done;

