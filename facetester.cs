using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using EyE.UnityAssetTypes;
using EyE.Threading;
using EyE.Geometry;

public class facetester : MonoBehaviour
{
    public FacesAndNeighbors faces;
    public MeshFilter mf;
    public GameObject facePRefab;
    private void OnEnable()
    {
        mf.sharedMesh = faces.meshRef;
        foreach (FaceDetails face in faces.faceDetails)
        {
            Instantiate(facePRefab, face.normal, Quaternion.FromToRotation(Vector3.forward, face.normal),transform);
        }
    }
}
