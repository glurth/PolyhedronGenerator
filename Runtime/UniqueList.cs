using System.Collections;
using System.Collections.Generic;
/// <summary>
/// Represents a unique list of elements. 
/// The list maintains uniqueness and provides O(1) lookup for Contains and IndexOf using an internal dictionary.
/// The cost is additional memory usage, and slower insert and remove operations.
/// Add and indexing operations remain O(1).
/// </summary>
public class UniqueList<K> : IList<K>
{
    List<K> internalList= new List<K>();                 // The underlying list to hold the elements, and keep them in order.
    Dictionary<K, int> internalDictionary= new Dictionary<K, int>(); // The dictionary for O(1) lookup by key.
    public UniqueList(){ }
    public UniqueList(IList<K> copyFrom)
    {
        foreach (K element in copyFrom)
        {
            Add(element);
        }
    }

    /// <summary>
    /// Rebuilds the internal dictionary from clear, by iterating through the list and assigning an index for each item.
    /// This method is used to ensure that the dictionary is up-to-date.
    /// </summary>
    void RebuildInternalDic()
    {
        internalDictionary.Clear();
        for (int i = 0; i < internalList.Count; i++)
            internalDictionary.Add(internalList[i], i);
    }
    
    /// <summary>
    /// Updates dictionary with correct index for all items in the list. 
    /// NOTE: ASSUMES that no entries in the dictionary exist that do not also exist in the internal list.  Such entires will not be updated.
    /// </summary>
    void UpdateInternalDicValues()
    {
        for (int i = 0; i < internalList.Count; i++)
            internalDictionary[internalList[i]]= i;
    }

    /// <summary>
    /// Indexer for accessing elements by index. O(1) for list access.
    /// The setter removes the old value from the dictionary and adds the new value.
    /// </summary>
    public K this[int index]
    {
        get => internalList[index]; // O(1) for list access.

        set
        {
            
            K oldValue = internalList[index]; // O(1) to access old value.
            if (internalDictionary.ContainsKey(value)) 
                throw new System.Exception("Cannot set elemnt in UniqueList to given value, it already exists in the list.");
            internalDictionary.Remove(oldValue); // O(1) for dictionary removal.
            internalList[index] = value; // O(1) for list modification.
            internalDictionary.Add(value, index); // O(1) for dictionary addition.
        }
    }

    /// <summary>
    /// Gets the count of elements in the list. O(1) as it's a property of the list.
    /// </summary>
    public int Count => internalList.Count;

    /// <summary>
    /// Indicates whether the list is read-only. O(1).
    /// </summary>
    public bool IsReadOnly => false;

    /// <summary>
    /// Adds an item to the list and updates the dictionary. 
    /// Throws an exception if the item already exists. 
    /// This operation is O(1) for both list and dictionary operations.
    /// </summary>
    public void Add(K item)
    {
        if (internalDictionary.ContainsKey(item))
            throw new System.ArgumentException("Item already exists in unique list");
        internalDictionary.Add(item, internalList.Count); // Add to dictionary as key with index it will have in the list as the value.
        internalList.Add(item);// Add to list.
    }

    /// <summary>
    /// Returns the index of the item in this list, after either being found-in or added-to the list.
    /// </summary>
    /// <param name="item">item to find or add</param>
    /// <returns></returns>
    public int GetIndexOrAdd(K item)
    {
        if (internalDictionary.TryGetValue(item, out int i))
            return i;
        internalDictionary.Add(item, internalList.Count); // Add to dictionary with new index.
        int idx = internalList.Count;
        internalList.Add(item);                           // Add to list.
        return idx;
    }


    /// <summary>
    /// Clears the list and the dictionary. O(n) for clearing list, O(1) for clearing dictionary.
    /// </summary>
    public void Clear()
    {
        internalList.Clear();
        internalDictionary.Clear();
    }

    /// <summary>
    /// Checks if an item exists in the list using the dictionary. 
    /// O(1) lookup.
    /// </summary>
    public bool Contains(K item)
    {
        return internalDictionary.ContainsKey(item);
    }

    /// <summary>
    /// Copies the elements to the specified array starting at the given index. O(n).
    /// </summary>
    public void CopyTo(K[] array, int arrayIndex)
    {
        internalList.CopyTo(array, arrayIndex);
    }

    /// <summary>
    /// Gets an enumerator for the list. O(1) for getting the enumerator.
    /// </summary>
    public IEnumerator<K> GetEnumerator()
    {
        return internalList.GetEnumerator();
    }

    /// <summary>
    /// Finds the index of an item using the dictionary. O(1) lookup.
    /// </summary>
    public int IndexOf(K item)
    {
        if (internalDictionary.TryGetValue(item, out int i))
            return i;
        return -1;
    }

    /// <summary>
    /// Inserts an item at the specified index.
    /// This operation is slow because it requires shifting elements in the list (O(n)) and rebuilding the dictionary (O(n)).
    /// </summary>
    public void Insert(int index, K item)
    {
        internalList.Insert(index, item); // O(n) for shifting elements.
        internalDictionary.Add(item, index);
        UpdateInternalDicValues();
        //RebuildInternalDic(); // O(n) to rebuild dictionary.
    }

    /// <summary>
    /// Removes an item from the list. 
    /// This operation is slow because it requires shifting elements in the list (O(n)) and rebuilding the dictionary (O(n)).
    /// </summary>
    public bool Remove(K item)
    {
        if (internalDictionary.TryGetValue(item, out int index))
        {
            internalList.RemoveAt(index);
            internalDictionary.Remove(item);
            UpdateInternalDicValues();
            return true;
        }
        return false;
    }

    /// <summary>
    /// Removes an item at the specified index. 
    /// This operation is slow because it requires shifting elements in the list (O(n)) and rebuilding the dictionary (O(n)).
    /// </summary>
    public void RemoveAt(int index)
    {
        if (index < internalList.Count && index >= 0)
        {
            K keyValue = internalList[index];
            internalList.RemoveAt(index);  // O(n) for shifting elements.
            internalDictionary.Remove(keyValue);
            UpdateInternalDicValues();
            //RebuildInternalDic(); // O(n) to rebuild dictionary.
        }
    }

    /// <summary>
    /// Gets an enumerator for the list, as required by IEnumerable. O(1) for getting the enumerator.
    /// </summary>
    IEnumerator IEnumerable.GetEnumerator()
    {
        return ((IEnumerable)internalList).GetEnumerator();
    }
}
