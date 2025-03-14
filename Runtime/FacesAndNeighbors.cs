using System.Collections;
using System.Collections.Generic;
using UnityEngine;
[System.Serializable]
public class FaceDetails
{
    public int index;
    [NonReorderable]
    public List<int> triangles = new List<int>();// index of triangles in this face, based upon the triangles list in the meshRef
    [NonReorderable]
    public List<int> neighborIndices = new List<int>();  //index of neighboring faces, based upon the list in FacesAndNeighbors
    public Vector3 normal;
}

public class FacesAndNeighbors : ScriptableObject
{
    public Mesh meshRef;
    [NonReorderable]
    public List<FaceDetails> faceDetails;

    public int GetFaceIndex(int meshTriangleIndex)
    {
        //search all faces
        for (int i = 0; i < faceDetails.Count; i++)
        {
            FaceDetails face = faceDetails[i];
            //check if in face's triangle list
            foreach (int triIndex in face.triangles)
            {
                if (triIndex == meshTriangleIndex)
                    return i;
            }
        }
        return -1;
    }

    //assume normal param, and normal in mesh have magnitude 1.0 (or, at least, all have the same magnitude)
    public int GetClosestFaceIndex(Vector3 normal)
    {
        int bestFace = -1;
        float closestDot = -1;
        for (int i = 0; i < faceDetails.Count; i++)
        {
            FaceDetails face = faceDetails[i];
            float howClose = Vector3.Dot(normal, face.normal);
            if (howClose > closestDot)
            {
                bestFace = i;
                closestDot = howClose;
            }
        }
        return bestFace;
    }
}
