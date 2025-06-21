using UnityEngine;
using EyE.Geometry;
namespace EyE.UnityAssetTypes
{

    public class PolyhedronSO : ScriptableObject
    {
        public Polyhedron poly;
        public FacesAndNeighbors facesAndNeighbors;
        public Mesh mesh;
    }
}