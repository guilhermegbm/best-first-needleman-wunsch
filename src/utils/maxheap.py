# Código original em https://www.geeksforgeeks.org/max-heap-in-python/
import sys

class MaxHeap:

    def __init__(self):
        self.heap = []
        self.FRONT = 0

    def size(self):
        return len(self.heap)

    def parent(self, pos):
        if pos == 0:
            return 0
        else:
            return (pos - 1) // 2

    def leftChild(self, pos):
        return (2 * pos) + 1

    def hasLeftChild(self, pos):
        return self.leftChild(pos) < self.size()

    def rightChild(self, pos):
        return (2 * pos) + 2

    def hasRightChild(self, pos):
        return self.rightChild(pos) < self.size()

    #def isLeaf(self, pos):
    #    if pos >= (self.size() // 2) and pos <= self.size():
    #        return True
    #    return False

    def swap(self, fpos, spos):
        self.heap[fpos], self.heap[spos] = (self.heap[spos], self.heap[fpos])

    def maxHeapify(self, pos):

        if ((self.hasLeftChild(pos) and self.heap[pos] < self.heap[self.leftChild(pos)])
            or (self.hasRightChild(pos) and self.heap[pos] < self.heap[self.rightChild(pos)])):

            if (self.hasRightChild(pos) and self.heap[self.rightChild(pos)] >= self.heap[self.leftChild(pos)]):
                self.swap(pos, self.rightChild(pos))
                self.maxHeapify(self.rightChild(pos))
            else:
                self.swap(pos, self.leftChild(pos))
                self.maxHeapify(self.leftChild(pos))

    def insert(self, element):
        self.heap.append(element)

        current = self.size() - 1

        while (self.heap[current] > self.heap[self.parent(current)]):
            self.swap(current, self.parent(current))
            current = self.parent(current)

    def isEmpty(self):
        return self.size() == 0

    def extractMax(self):

        if self.isEmpty():
            raise Exception("Empty Heap")

        first = self.heap[self.FRONT]
        last = self.heap.pop()

        if self.size() > 0:
            self.heap[self.FRONT] = last
            self.maxHeapify(self.FRONT)

        return first

    def Print(self):
        print(self.heap)