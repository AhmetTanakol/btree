#ifndef INCLUDE_MODERNDBS_BTREE_H
#define INCLUDE_MODERNDBS_BTREE_H

#include <cstddef>
#include <cstring>
#include <functional>
#include <experimental/optional>
#include <iostream>
#include <vector>
#include "moderndbs/buffer_manager.h"
#include "moderndbs/defer.h"
#include "moderndbs/segment.h"
#include "moderndbs/defer.h"

using Defer = moderndbs::Defer;

namespace moderndbs {

    template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
    struct BTree : public Segment {
        struct Node {
            /// The level in the tree.
            uint16_t level;

            /// The number of children.
            uint16_t count;

            // Constructor
            Node(uint16_t level, uint16_t count)
                    : level(level), count(count) {}

            /// Is the node a leaf node?
            bool is_leaf() const { return level == 0; }

            std::experimental::optional<uint64_t> parentPageId;
        };

        struct InnerNode: public Node {
            /// The capacity of a node.
            /// TODO think about the capacity that the nodes have.
            static constexpr uint32_t kCapacity = (PageSize / (sizeof(KeyT) + sizeof(ValueT))) - 2;

            /// The keys.
            KeyT keys[kCapacity];

            /// The children.
            uint64_t children[kCapacity + 1];

            /// Constructor.
            InnerNode() : Node(0, 0) {}

            /// Binary search implementation for lower bound
            /// Get the index of the first key that is not less than than a provided key.
            /// @param[in] key          The key that should be searched.
            std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
                std::pair<uint32_t , bool> result;
                int low = 0;
                int high = this->count - 1;
                std::experimental::optional<uint32_t> index;
                if (high == 1 && !ComparatorT()(this->keys[0], key)) {
                    result.first = 0;
                    result.second = true;
                    return result;
                }
                while (low <= high) {
                    int mid = low + (high -low) / 2;

                    if (key == this->keys[mid]) {
                        result.first = static_cast<uint32_t>(mid);
                        result.second = true;
                        return result;
                    } else if (key > this->keys[mid]) {
                        low = mid + 1;
                    } else {
                        index = mid;
                        high = mid - 1;
                    }
                }
                if (index) {
                    result.first = *index;
                    result.second = true;
                    return result;
                }
                result.first = 0;
                result.second = false;
                return result;
            }
            /// Get the index of the first key that is not less than than a provided key.
            /// @param[in] key          The key that should be searched.
            /*std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
                std::pair<uint32_t , bool> result;
                for(int i=0;i < this->count;i++) {
                    if(!ComparatorT()(this->keys[i], key)) {
                        result.first = i;
                        result.second = true;
                        return result;
                    }
                }
                result.first = 0;
                result.second = false;
                return result;
            }*/

            /// Insert a key.
            /// @param[in] key          The separator that should be inserted.
            /// @param[in] split_page   The id of the split page that should be inserted.
            void insert(const KeyT &key, uint64_t split_page, const std::experimental::optional<uint64_t> &leftNode) {
                std::vector<KeyT> tempKeys;
                std::vector<uint64_t> tempChildren;
                std::experimental::optional<int> index;
                for (int i=0;i < this->count - 1;i++) {
                    tempKeys.push_back(this->keys[i]);
                    tempChildren.push_back(this->children[i]);
                    if(!ComparatorT()(this->keys[i], key) && !index) {
                        index = i;
                    }
                    this->keys[i] = 0;
                    this->children[i] = 0;
                }
                /// push the last element in the children array
                tempChildren.push_back(this->children[this->count - 1]);
                this->children[this->count - 1] = 0;
                if (!index) {
                    tempKeys.insert(tempKeys.begin() + (this->count - 1), key);
                    tempChildren.insert(tempChildren.begin() + this->count, split_page);
                } else {
                    tempKeys.insert(tempKeys.begin() + *index, key);
                    if (leftNode) {
                        index = *index + 1;
                    }
                    tempChildren.insert(tempChildren.begin() + *index, split_page);
                }
                for (int i=0;i < static_cast<int>(tempKeys.size()); i++) {
                    this->keys[i] = tempKeys[i];
                }
                for (int i=0;i < static_cast<int>(tempChildren.size()); i++) {
                    this->children[i] = tempChildren[i];
                }
                this->count++;
            }

            /// Split the node.
            /// @param[in] buffer       The buffer for the new page.
            /// @return                 The separator key.
            KeyT split(std::byte* buffer) {
                auto mid_index = 1 + ((this->kCapacity - 1) / 2);
                auto new_inner_node = reinterpret_cast<InnerNode*>(buffer);
                int startIndex = 0;
                /// calculate a separator key
                auto separator_key = static_cast<int>((this->keys[mid_index] + this->keys[mid_index - 1]) / 2);
                for(int i=mid_index;i < static_cast<int>(this->kCapacity) + 1;i++) {
                    if (startIndex == 0) {
                        new_inner_node->keys[startIndex] = this->keys[i];
                        new_inner_node->children[startIndex] = this->children[i];
                        this->keys[i] = 0;
                        this->children[i] = 0;
                        this->count--;
                        startIndex++;
                        i++;
                        new_inner_node->children[startIndex] = this->children[i];
                        this->children[i] = 0;
                        this->count--;
                        startIndex++;
                    } else {
                        new_inner_node->keys[startIndex - 1] = this->keys[i - 1];
                        new_inner_node->children[startIndex] = this->children[i];
                        this->keys[i - 1] = 0;
                        this->children[i] = 0;
                        this->count--;
                        startIndex++;
                    }
                }
                new_inner_node->count = startIndex;
                return separator_key;
            }
        };

        struct LeafNode: public Node {
            /// The capacity of a node.
            /// TODO think about the capacity that the nodes have.
            static constexpr uint32_t kCapacity = (PageSize / (sizeof(KeyT) + sizeof(ValueT))) - 2;

            /// The keys.
            KeyT keys[kCapacity];

            /// The values.
            ValueT values[kCapacity];

            /// Constructor.
            LeafNode() : Node(0, 0) {}

            /// Insert a key.
            /// @param[in] key          The key that should be inserted.
            /// @param[in] value        The value that should be inserted.
            void insert(const KeyT &key, const ValueT &value) {
                if (this->count == 0) {
                    this->keys[this->count] = key;
                    this->values[this->count] = value;
                    this->count++;
                } else {
                    /// sort keys and shift keys (including values) wrt incoming key
                    std::vector<KeyT> tempKeys;
                    std::vector<KeyT> tempValues;
                    std::experimental::optional<int> index;
                    for (int i=0;i < this->count;i++) {
                        tempKeys.push_back(this->keys[i]);
                        tempValues.push_back(this->values[i]);
                        if(!ComparatorT()(this->keys[i], key) && !index) {
                            index = i;
                        }
                        this->keys[i] = 0;
                        this->values[i] = 0;
                    }
                    if (!index) {
                        index = this->count;
                        /// if the key already exists, overwrites its value
                        if (tempKeys[*index - 1] == key) {
                            tempValues[*index - 1] = value;
                        } else {
                            tempKeys.insert(tempKeys.begin() + *index, key);
                            tempValues.insert(tempValues.begin() + *index, value);
                            this->count++;
                        }
                    } else {
                        if (tempKeys[*index] == key) {
                            tempValues[*index] = value;
                        } else {
                            tempKeys.insert(tempKeys.begin() + *index, key);
                            tempValues.insert(tempValues.begin() + *index, value);
                            this->count++;
                        }
                    }
                    for (int i=0;i < static_cast<int>(tempKeys.size()); i++) {
                        this->keys[i] = tempKeys[i];
                        this->values[i] = tempValues[i];
                    }
                }
            }

            /// Erase a key.
            void erase(int &index) {
                // delete add shift other elements to left
                for (int i = index; i < this->count; ++i) {
                    this->keys[i] = this->keys[i + 1];
                    this->values[i] = this->values[i + 1];
                }
                this->keys[this->count] = 0;
                this->values[this->count] = 0;
                this->count--;
            }

            /// Split the node.
            /// @param[in] buffer       The buffer for the new page.
            /// @return                 The separator key.
            KeyT split(std::byte* buffer) {
                auto mid_index = 1 + ((this->kCapacity - 1) / 2);
                auto new_leaf_node = reinterpret_cast<LeafNode*>(buffer);
                int startIndex = 0;
                for(int i=mid_index;i < static_cast<int>(this->kCapacity);i++) {
                    new_leaf_node->keys[startIndex] = this->keys[i];
                    new_leaf_node->values[startIndex] = this->values[i];
                    this->keys[i] = 0;
                    this->values[i] = 0;
                    this->count--;
                    startIndex++;
                }
                new_leaf_node->count = startIndex;
                return this->keys[startIndex - 1];
            }
        };

        /// The root.
        std::experimental::optional<uint64_t> root;

        bool isTreeEmpty;

        uint16_t rootLevel = 0;

        std::map<KeyT, bool> deletedKeys;

        /// Next page id.
        /// You don't need to worry about about the page allocation.
        /// (Neither fragmentation, nor persisting free-space bitmaps)
        /// Just increment the next_page_id whenever you need a new page.
        uint64_t next_page_id;

        /// Constructor.
        BTree(uint16_t segment_id, BufferManager &buffer_manager)
                : Segment(segment_id, buffer_manager) {
            this->isTreeEmpty = true;
        }

        /// Lookup an entry in the tree.
        /// @param[in] key      The key that should be searched.
        std::experimental::optional<ValueT> lookup(const KeyT &key) {
            bool keyFound = false;
            std::experimental::optional<ValueT> foundKey;
            if (this->deletedKeys.find(key) != this->deletedKeys.end()) {
                return foundKey;
            }
            if (!this->root) {
                return foundKey;
            }
            auto current_page_id = *this->root;
            int next = 0;
            uint64_t previousParentPageId = current_page_id;
            while (!keyFound) {
                auto& current_page = this->buffer_manager.fix_page(current_page_id, false);
                auto current_node = reinterpret_cast<Node*>(current_page.get_data());
                if (current_node->is_leaf()) {
                    auto leaf_node = reinterpret_cast<LeafNode*>(current_node);
                    for (int i=0;i < leaf_node->count; i++) {
                        if (leaf_node->keys[i] == key) {
                            foundKey = leaf_node->values[i];
                            return foundKey;
                        }
                    }
                    /// In this implementation, we can find any key at most 2 searches
                    /// at this point, no need to check further
                    if (next == 2) {
                        return foundKey;
                    }
                    /// we need to go back to parent node and check the next children of it
                    current_page_id = previousParentPageId;
                    next++;
                } else {
                    auto inner_node = reinterpret_cast<InnerNode*>(current_node);
                    if (current_node->parentPageId) {
                        previousParentPageId = *current_node->parentPageId;
                    }
                    /// move to next node
                    std::pair<uint32_t, bool> lower_bound = inner_node->lower_bound(key);
                    if (lower_bound.second) {
                        if (inner_node->count < next) {
                            return foundKey;
                        } else {
                            current_page_id = inner_node->children[lower_bound.first + next];
                            next = 0;
                        }
                    } else {
                        current_page_id = inner_node->children[inner_node->count - 1];
                    }
                }
                this->buffer_manager.unfix_page(current_page, false);
            }
            return foundKey;
        }

        /// Erase an entry in the tree.
        /// @param[in] key      The key that should be searched.
        void erase(const KeyT &key) {
            if (this->root) {
                bool keyFound = false;
                auto current_page_id = *this->root;
                int next = 0;
                uint64_t previousParentPageId = current_page_id;
                while (!keyFound) {
                    auto& current_page = this->buffer_manager.fix_page(current_page_id, true);
                    auto current_node = reinterpret_cast<Node*>(current_page.get_data());
                    if (current_node->is_leaf()) {
                        auto leaf_node = reinterpret_cast<LeafNode*>(current_node);
                        for (int i=0;i < leaf_node->count; i++) {
                            if (leaf_node->keys[i] == key) {
                                leaf_node->erase(i);
                                this->deletedKeys.insert(std::pair<KeyT,bool>(key,true));
                                keyFound = true;
                                break;
                            }
                        }
                        current_page_id = previousParentPageId;
                        next++;
                    } else {
                        auto inner_node = reinterpret_cast<InnerNode*>(current_node);
                        if (current_node->parentPageId) {
                            previousParentPageId = *current_node->parentPageId;
                        }
                        /// move to next node
                        std::pair<uint32_t, bool> lower_bound = inner_node->lower_bound(key);
                        if (lower_bound.second) {
                            if (inner_node->count < next) {
                                keyFound = true;
                            } else {
                                current_page_id = inner_node->children[lower_bound.first + next];
                                next = 0;
                            }
                        } else {
                            current_page_id = inner_node->children[inner_node->count - 1];
                        }
                    }
                    this->buffer_manager.unfix_page(current_page, true);
                }
            }
        }

        /// Inserts a new entry into the tree.
        /// @param[in] key      The key that should be inserted.
        /// @param[in] value    The value that should be inserted.
        void insert(const KeyT &key, const ValueT &value) {
            if (this->isTreeEmpty) {
                this->root = 0;
                this->next_page_id = 1;
                this->isTreeEmpty = false;
            }
            /// use this variable to iterate over multiple inner nodes
            auto current_node_page_id = *this->root;
            bool isKeyInserted = false;
            //std::cout << "key to insert: " << key << '\n';
            while (!isKeyInserted) {
                auto& current_page = this->buffer_manager.fix_page(current_node_page_id, true);
                auto current_node = reinterpret_cast<Node*>(current_page.get_data());
                /// if root node is a leaf node
                if (current_node->is_leaf()) {
                    auto leaf_node = static_cast<LeafNode*>(current_node);
                    /// check if there is space for the key
                    if (leaf_node->count < leaf_node->kCapacity) {
                        leaf_node->insert(key, value);
                        this->buffer_manager.unfix_page(current_page, true);
                        isKeyInserted = true;
                    } else {
                        /// create new leaf node
                        auto new_leaf_page_id = this->next_page_id;
                        this->next_page_id++;
                        auto& new_leaf_page = this->buffer_manager.fix_page(new_leaf_page_id, true);
                        KeyT separator_key = leaf_node->split(reinterpret_cast<std::byte *>(new_leaf_page.get_data()));
                        auto new_node = reinterpret_cast<Node*>(new_leaf_page.get_data());
                        auto new_leaf_node = static_cast<LeafNode*>(new_node);
                        /// add new key
                        /// if separator_key is not lower than key
                        /// add to lhs, otherwise add to rhs
                        if(!ComparatorT()(separator_key, key)) {
                            leaf_node->insert(key, value);
                        } else {
                            new_leaf_node->insert(key, value);
                        }
                        /*std::cout << "leaf node: " << current_node_page_id << '\n';
                        for (int i=0;i < leaf_node->count;i++) {
                            std::cout << leaf_node->keys[i] << ' ';
                        }
                        std::cout << '\n';
                        std::cout << "new leaf node: " << new_leaf_page_id << '\n';
                        for (int i=0;i < new_leaf_node->count;i++) {
                            std::cout << new_leaf_node->keys[i] << ' ';
                        }
                        std::cout << '\n';*/

                        /// check if current node has a parent or not
                        if (!leaf_node->parentPageId) {
                            /// root has new page id
                            this->root = this->next_page_id;
                            this->next_page_id++;
                            auto& new_root_node_page = this->buffer_manager.fix_page(*this->root, true);
                            auto new_inner_root_node = reinterpret_cast<Node *>(new_root_node_page.get_data());
                            auto new_root_node = static_cast<InnerNode *>(new_inner_root_node);
                            /// increase level by 1
                            rootLevel += 1;
                            new_root_node->level = rootLevel;
                            /// add first key
                            new_root_node->keys[0] = separator_key;
                            /// add children to first key
                            new_root_node->children[0] = current_node_page_id;
                            new_root_node->children[1] = new_leaf_page_id;
                            /// increase children number
                            new_root_node->count += 2;
                            /// set parent page id for both leaf nodes
                            leaf_node->parentPageId = *this->root;
                            new_leaf_node->parentPageId = *this->root;
                            /*std::cout << "creating new parent for leaf node, page id: " << *this->root << std::endl;
                            std::cout << "root node" << '\n';
                            for (int i=0;i < new_root_node->count - 1;i++) {
                                std::cout << new_root_node->keys[i] << ' ';
                            }
                            std::cout << '\n';
                            std::cout << "----------------------------------------------"  << std::endl;*/

                            this->buffer_manager.unfix_page(new_leaf_page, true);
                            this->buffer_manager.unfix_page(new_root_node_page, true);
                        } else {
                            /// add separator_key into existing parent node
                            auto& parent_node_page = this->buffer_manager.fix_page(*leaf_node->parentPageId, true);
                            auto parent_node = reinterpret_cast<Node *>(parent_node_page.get_data());
                            auto parent_inner_node = static_cast<InnerNode *>(parent_node);
                            //std::cout << " adding into existing parent node for leaf, page id: " << *leaf_node->parentPageId << " root level: " << parent_inner_node->level << ", " << " children count: " << parent_inner_node->count << std::endl;
                            std::experimental::optional<uint64_t> leftNode;
                            if (leaf_node->values[leaf_node->count -1] < new_leaf_node->values[0]) {
                                leftNode = current_node_page_id;
                            }
                            parent_inner_node->insert(separator_key, new_leaf_page_id, leftNode);
                            /*std::cout << "root node" << '\n';
                            std::cout << "root node keys" << '\n';
                            int numberOfKeys = 1;
                            int numberOfChild = 1;
                            for (int i=0;i < parent_inner_node->count - 1;i++) {
                                std::cout << parent_inner_node->keys[i] << ' ';
                                numberOfKeys++;
                            }
                            std::cout << '\n';
                            std::cout << "root node children" << '\n';
                            for (int i=0;i < parent_inner_node->count;i++) {
                                std::cout << parent_inner_node->children[i] << ' ';
                                numberOfChild++;
                            }
                            std::cout << '\n';
                            std::cout << "#keys: " << numberOfKeys << " #children: " << numberOfChild << '\n';
                            if (numberOfKeys == 64) {
                                std::cout << "asd" << std::endl;
                            }
                            std::cout << '\n';
                            std::cout << "*************************************************************************" << '\n';*/
                            new_leaf_node->parentPageId = *leaf_node->parentPageId;
                            this->buffer_manager.unfix_page(parent_node_page, true);
                        }
                        this->buffer_manager.unfix_page(current_page, true);
                        isKeyInserted = true;
                    }
                } else {
                    auto inner_node = static_cast<InnerNode*>(current_node);
                    /// create new inner node if current inner node is full
                    if (inner_node->count == (inner_node->kCapacity + 1)) {
                        /// create new inner node
                        auto new_inner_node_page_id = this->next_page_id;
                        this->next_page_id++;
                        auto& new_inner_node_page = this->buffer_manager.fix_page(new_inner_node_page_id, true);
                        KeyT separator_key = inner_node->split(reinterpret_cast<std::byte *>(new_inner_node_page.get_data()));
                        auto new_node = reinterpret_cast<Node*>(new_inner_node_page.get_data());
                        auto new_inner_node = static_cast<InnerNode*>(new_node);
                        new_inner_node->level = inner_node->level;
                        /// set the new parent id of the children of the new inner node
                        for (int i=0;i < new_inner_node->count;i++) {
                            auto& child = this->buffer_manager.fix_page(new_inner_node->children[i], true);
                            auto child_node = reinterpret_cast<Node*>(child.get_data());
                            auto child_inner_node = static_cast<InnerNode*>(child_node);
                            child_inner_node->parentPageId = new_inner_node_page_id;
                        }
                        /// check if current inner node has a parent
                        if (!inner_node->parentPageId) {
                            /// root has new page id
                            this->root = this->next_page_id;
                            this->next_page_id++;
                            auto& new_root_node_page = this->buffer_manager.fix_page(*this->root, true);
                            auto new_root_node = reinterpret_cast<Node *>(new_root_node_page.get_data());
                            auto new_root_inner_node = static_cast<InnerNode *>(new_root_node);
                            /// increase level by 1
                            rootLevel += 1;
                            new_root_inner_node->level = rootLevel;
                            /// add first key
                            new_root_inner_node->keys[0] = separator_key;
                            /// add children to first key
                            new_root_inner_node->children[0] = current_node_page_id;
                            new_root_inner_node->children[1] = new_inner_node_page_id;
                            /// increase children number
                            new_root_inner_node->count += 2;
                            /// set parent page id for both inner nodes
                            inner_node->parentPageId = *this->root;
                            new_inner_node->parentPageId = *this->root;
                            /*std::cout << "creating new parent for inner node, page id: " << *this->root << std::endl;
                            std::cout << "inner node" << '\n';
                            for (int i=0;i < inner_node->count;i++) {
                                std::cout << inner_node->keys[i] << ' ';
                            }
                            std::cout << '\n';
                            std::cout << "new inner node" << '\n';
                            for (int i=0;i < new_inner_node->count - 1;i++) {
                                std::cout << new_inner_node->keys[i] << ' ';
                            }
                            std::cout << '\n';
                            std::cout << "root node" << '\n';
                            for (int i=0;i < new_root_inner_node->count - 1;i++) {
                                std::cout << new_root_inner_node->keys[i] << ' ';
                            }
                            std::cout << '\n';
                            std::cout << "----------------------------------------------"  << std::endl;*/

                            this->buffer_manager.unfix_page(new_root_node_page, true);
                            /// move to next node
                            std::pair<uint32_t, bool> lower_bound = new_root_inner_node->lower_bound(key);
                            if (lower_bound.second) {
                                /// go to lhs child
                                current_node_page_id = new_root_inner_node->children[lower_bound.first];
                            } else {
                                /// go to rhs child
                                current_node_page_id = new_root_inner_node->children[new_root_inner_node->count - 1];
                            }
                        } else {
                            /// get the parent node of the current node
                            auto& parent_node_page = this->buffer_manager.fix_page(*inner_node->parentPageId, true);
                            auto parent_node = reinterpret_cast<Node *>(parent_node_page.get_data());
                            auto parent_inner_node = static_cast<InnerNode *>(parent_node);
                            //std::cout << "current page id: " << current_node_page_id << " adding into existing parent for inner, page id: " << *inner_node->parentPageId << " root level: " << parent_inner_node->level << " children count: " << parent_inner_node->count << std::endl;
                            std::experimental::optional<uint64_t> leftNode;
                            parent_inner_node->insert(separator_key, new_inner_node_page_id, leftNode);
                            new_inner_node->parentPageId = *inner_node->parentPageId;
                            /*std::cout << "root node" << '\n';
                            for (int i=0;i < parent_inner_node->count;i++) {
                                std::cout << parent_inner_node->keys[i] << ' ';
                            }
                            std::cout << '\n';*/
                            /// move to next node
                            std::pair<uint32_t, bool> lower_bound = parent_inner_node->lower_bound(key);
                            if (lower_bound.second) {
                                /// go to lhs child
                                current_node_page_id = parent_inner_node->children[lower_bound.first];
                            } else {
                                /// go to rhs child
                                current_node_page_id = parent_inner_node->children[parent_inner_node->count - 1];
                            }
                            this->buffer_manager.unfix_page(parent_node_page, true);
                        }
                        this->buffer_manager.unfix_page(new_inner_node_page, true);
                    } else {
                        /// move to next node
                        /// if key is greater than any of the keys in the inner node
                        /// go to last children
                        std::pair<uint32_t, bool> lower_bound = inner_node->lower_bound(key);
                        if (lower_bound.second) {
                            /// go to lhs child
                            current_node_page_id = inner_node->children[lower_bound.first];
                        } else {
                            /// go to rhs child
                            current_node_page_id = inner_node->children[inner_node->count - 1];
                        }
                    }
                    this->buffer_manager.unfix_page(current_page, true);
                }
            }
        }
    };

}  // namespace moderndbs

#endif
