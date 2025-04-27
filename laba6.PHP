<?php
function customSort(array $array): array {
    $evenIndexes = [];
    $evenValues = [];
    $oddIndexes = [];
    $oddValues = [];

    foreach ($array as $i => $value) {
        if ($value === 0) continue;
        
        if ($value % 2 === 0) {
            $evenIndexes[] = $i;
            $evenValues[] = $value;
        } else {
            $oddIndexes[] = $i;
            $oddValues[] = $value;
        }
    }

    usort($oddValues, function($a, $b) {
        return $a <=> $b; 
    });
    
    usort($evenValues, function($a, $b) {
        return $b <=> $a; 
    });

    $sorted = $array;
    
    foreach ($oddIndexes as $k => $index) {
        $sorted[$index] = $oddValues[$k];
    }
    
    foreach ($evenIndexes as $k => $index) {
        $sorted[$index] = $evenValues[$k];
    }

    return $sorted;
}

$array = [5, 1, 2, 9, 0, 8, 7, 0, 6];
$result = customSort($array);
echo "Результат: " . implode(', ', $result); 
?>
